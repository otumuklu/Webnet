#!/usr/bin/env python3
"""Refresh the public aerospace briefing from NASA, SpaceNews, arXiv and Grants.gov."""

from __future__ import annotations

import json
import re
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta, timezone
from email.utils import parsedate_to_datetime
from html.parser import HTMLParser
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode, urlparse
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
SOURCE_OUTPUT = ROOT / "links" / "auto-news.json"
SITE_OUTPUT = ROOT / "docs" / "links" / "auto-news.json"
USER_AGENT = "HAVA-Lab-News/1.0 (mailto:tumuko@rpi.edu)"
INDUSTRY_ITEMS_PER_SOURCE = 2
PAPER_SOURCE_LIMITS = {
    "arXiv": 2,
    "AIAA": 3,
    "Physics of Fluids": 2,
    "Physical Review Fluids": 2,
    "Journal of Fluid Mechanics": 3,
}


class TextExtractor(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.parts: list[str] = []
        self.hidden_depth = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag in {"script", "style"}:
            self.hidden_depth += 1

    def handle_endtag(self, tag: str) -> None:
        if tag in {"script", "style"} and self.hidden_depth:
            self.hidden_depth -= 1

    def handle_data(self, data: str) -> None:
        if not self.hidden_depth:
            self.parts.append(data)


class BlueOriginCardParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.cards: list[dict[str, str]] = []
        self.card: dict[str, object] | None = None
        self.anchor_depth = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        if tag == "li" and self.card is None:
            self.card = {"text": [], "links": []}
        if self.card is not None and tag == "a":
            href = attributes.get("href", "") or ""
            if href.startswith("/news/"):
                self.card["links"].append({"url": href, "title": []})
                self.anchor_depth = len(self.card["links"])

    def handle_endtag(self, tag: str) -> None:
        if tag == "a":
            self.anchor_depth = 0
        if tag == "li" and self.card is not None:
            self.cards.append(self.card)
            self.card = None
            self.anchor_depth = 0

    def handle_data(self, data: str) -> None:
        if self.card is None:
            return
        self.card["text"].append(data)
        if self.anchor_depth:
            self.card["links"][self.anchor_depth - 1]["title"].append(data)


class SbirTopicListParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.topics: list[tuple[str, str]] = []
        self.current: tuple[str, list[str]] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "a":
            href = dict(attrs).get("href", "") or ""
            if re.fullmatch(r"/topics/\d+", href):
                self.current = (href, [])

    def handle_data(self, data: str) -> None:
        if self.current is not None:
            self.current[1].append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag == "a" and self.current is not None:
            href, title_parts = self.current
            self.topics.append((href, plain_text(" ".join(title_parts))))
            self.current = None


def plain_text(value: str) -> str:
    extractor = TextExtractor()
    extractor.feed(value)
    return re.sub(r"\s+", " ", " ".join(extractor.parts)).strip()


def safe_url(value: str) -> str:
    parsed = urlparse(value.strip())
    return value.strip() if parsed.scheme == "https" and parsed.netloc else ""


def fetch(url: str, headers: dict[str, str] | None = None) -> bytes:
    request_headers = {"User-Agent": USER_AGENT}
    if headers:
        request_headers.update(headers)
    request = Request(url, headers=request_headers)
    with urlopen(request, timeout=25) as response:
        return response.read()


def child_text(item: ET.Element, name: str) -> str:
    element = item.find(name)
    return "" if element is None or element.text is None else element.text.strip()


def rss_items(url: str, source: str) -> list[dict[str, str]]:
    root = ET.fromstring(fetch(url))
    items = []
    for item in root.findall(".//item"):
        title = plain_text(child_text(item, "title"))
        link = safe_url(child_text(item, "link"))
        description = plain_text(child_text(item, "description"))
        description = re.sub(r"\s*The post\s+.+?\s+appeared first on\s+.+\.?$", "", description, flags=re.IGNORECASE)
        date_text = child_text(item, "pubDate")
        publisher = plain_text(child_text(item, "source")) if "news.google.com" in url else ""
        if not title or not link:
            continue
        if "news.google.com" in url:
            description = description.replace(title, "", 1).replace(publisher, "", 1).strip(" -|:")
        try:
            date = parsedate_to_datetime(date_text).date().isoformat()
        except (TypeError, ValueError, OverflowError):
            date = ""
        items.append({"title": title, "url": link, "summary": description, "date": date, "source": publisher or source})
    return items


def paper_items() -> list[dict[str, str]]:
    papers = []
    for category in ("physics.flu-dyn", "physics.ao-ph"):
        root = ET.fromstring(fetch(f"https://rss.arxiv.org/rss/{category}"))
        for item in root.findall(".//item"):
            title = plain_text(child_text(item, "title"))
            link = safe_url(child_text(item, "link"))
            description = plain_text(child_text(item, "description"))
            date_text = child_text(item, "pubDate")
            try:
                date = parsedate_to_datetime(date_text).date().isoformat()
            except (TypeError, ValueError, OverflowError):
                date = ""
            summary = re.sub(r"^.*?Abstract:\s*", "", description, flags=re.IGNORECASE)
            summary = " ".join(re.split(r"(?<=[.!?])\s+", summary)[:2])
            if title and link:
                papers.append({"title": title, "url": link, "summary": summary[:360], "date": date, "source": "arXiv"})
    papers = deduplicate(papers)[:PAPER_SOURCE_LIMITS["arXiv"]]
    papers.extend(crossref_items("member:1387", "AIAA"))
    papers.extend(crossref_items("issn:1070-6631", "Physics of Fluids"))
    papers.extend(crossref_items("issn:2469-990X", "Physical Review Fluids"))
    papers.extend(crossref_items("issn:0022-1120", "Journal of Fluid Mechanics"))
    papers.sort(key=lambda item: item["date"], reverse=True)
    return deduplicate(papers)[:sum(PAPER_SOURCE_LIMITS.values())]


def crossref_items(journal_filter: str, source: str) -> list[dict[str, str]]:
    start_date = (datetime.now(timezone.utc).date() - timedelta(days=45)).isoformat()
    query = urlencode(
        {
            "filter": f"{journal_filter},from-pub-date:{start_date}",
            "rows": 100 if source == "AIAA" else 30,
            "select": "title,URL,published,published-online,abstract,container-title,publisher",
        }
    )
    url = f"https://api.crossref.org/works?{query}"
    try:
        for attempt in range(3):
            try:
                response = json.loads(fetch(url, {"Accept": "application/json"}))
                break
            except HTTPError as error:
                if error.code != 429 or attempt == 2:
                    raise
                retry_after = error.headers.get("Retry-After", "")
                delay = int(retry_after) if retry_after.isdigit() else 2 ** (attempt + 1)
                time.sleep(min(delay, 10))
    except (URLError, TimeoutError, OSError, json.JSONDecodeError) as error:
        print(f"Could not refresh {source} papers: {error}", file=sys.stderr)
        return []

    allowed_aiaa_titles = {
        "aiaa journal",
        "journal of aircraft",
        "journal of guidance, control, and dynamics",
        "journal of propulsion and power",
        "journal of spacecraft and rockets",
        "journal of thermophysics and heat transfer",
        "journal of aerospace information systems",
    }
    results = []
    today = datetime.now(timezone.utc).date().isoformat()
    for work in response.get("message", {}).get("items", []):
        container = (work.get("container-title") or [""])[0]
        if source == "AIAA" and container.lower() not in allowed_aiaa_titles:
            continue
        title = plain_text((work.get("title") or [""])[0])
        link = safe_url(str(work.get("URL", "")))
        published = work.get("published-online") or work.get("published") or {}
        date_parts = (published.get("date-parts") or [[]])[0]
        if len(date_parts) < 2:
            continue
        date = f"{date_parts[0]:04d}-{date_parts[1]:02d}"
        if len(date_parts) > 2:
            date += f"-{date_parts[2]:02d}"
        if date > today:
            continue
        abstract = plain_text(str(work.get("abstract", "")))
        summary = " ".join(re.split(r"(?<=[.!?])\s+", abstract)[:2])
        label = container if source == "AIAA" else source
        if title and link:
            results.append({"title": title, "url": link, "summary": summary[:420], "date": date, "source": label})
    results.sort(key=lambda item: item["date"], reverse=True)
    return results[:PAPER_SOURCE_LIMITS[source]]


def spacex_items() -> list[dict[str, str]]:
    updates = json.loads(fetch("https://content.spacex.com/api/spacex-website/updates"))
    items = []
    for update in updates:
        title = plain_text(str(update.get("title", "")))
        update_id = str(update.get("updateId", ""))
        date = str(update.get("date", ""))
        paragraphs = [
            plain_text(str(block.get("paragraph", "")))
            for block in update.get("contentBlocks", [])
            if block.get("paragraph")
        ]
        if title and update_id:
            items.append({
                "title": title,
                "url": f"https://www.spacex.com/updates#{quote(update_id)}",
                "summary": " ".join(paragraphs[:2])[:420],
                "date": date,
                "source": "SpaceX",
            })
    items.sort(key=lambda item: item["date"], reverse=True)
    return items[:INDUSTRY_ITEMS_PER_SOURCE]


def blue_origin_items() -> list[dict[str, str]]:
    parser = BlueOriginCardParser()
    html = fetch(
        "https://www.blueorigin.com/news",
        {
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/130.0.0.0 Safari/537.36",
            "Accept": "text/html,application/xhtml+xml",
            "Accept-Language": "en-US,en;q=0.9",
        },
    ).decode("utf-8", "replace")
    # Parse the official article cards while retaining their publisher-provided dates and summaries.
    parser.feed(html)
    items = []
    for card in parser.cards:
        card_text = plain_text(" ".join(card["text"]))
        date_match = re.search(r"\b[A-Z][a-z]{2}\s+\d{1,2},\s+\d{4}\b", card_text)
        try:
            date = datetime.strptime(date_match.group(), "%b %d, %Y").date().isoformat() if date_match else ""
        except ValueError:
            date = ""
        for link in card["links"]:
            title = plain_text("".join(link["title"]))
            href = safe_url("https://www.blueorigin.com" + link["url"])
            if not title or not href:
                continue
            summary = card_text.replace(date_match.group(), "", 1) if date_match else card_text
            summary = summary.replace(title, "", 1).replace("News", "", 1).strip()
            items.append({"title": title, "url": href, "summary": summary[:420], "date": date, "source": "Blue Origin"})
    items.sort(key=lambda item: item["date"], reverse=True)
    return deduplicate(items)[:INDUSTRY_ITEMS_PER_SOURCE]


def major_company_items() -> list[dict[str, str]]:
    searches = (
        (
            '("Lockheed Martin" OR Raytheon OR RTX OR "Northrop Grumman" OR '
            'L3Harris OR "Sierra Nevada" OR "Sierra Space" OR "Rocket Lab" OR '
            '"Firefly Aerospace" OR "Intuitive Machines" OR "Axiom Space") '
            '(aerospace OR space)',
            (
                ("Lockheed Martin",),
                ("Raytheon", "RTX"),
                ("Northrop Grumman",),
                ("L3Harris",),
                ("Sierra Nevada", "Sierra Space"),
                ("Rocket Lab",),
                ("Firefly Aerospace",),
                ("Intuitive Machines",),
                ("Axiom Space",),
            ),
        ),
        (
            '(Boeing OR Airbus) (aircraft OR aerospace)',
            (("Boeing",), ("Airbus",)),
        ),
        (
            '(JAXA OR "Japan Aerospace Exploration Agency" OR ESA OR '
            '"European Space Agency") space',
            (("JAXA", "Japan Aerospace Exploration Agency"), ("ESA", "European Space Agency")),
        ),
    )
    selected = []
    seen_urls: set[str] = set()
    excluded_sources = re.compile(
        r"britannica|encyclopedia|motley fool|yahoo finance|investorplace|"
        r"tradingview|marketwatch|stock titan|morningstar|investing\.com|finance|"
        r"kalkine|dars\.gov\.et|guru.?focus|ad hoc news|stock|ticker",
        re.IGNORECASE,
    )
    for query, organizations in searches:
        url = "https://news.google.com/rss/search?" + urlencode(
            {"q": query, "hl": "en-US", "gl": "US", "ceid": "US:en"}
        )
        candidates = rss_items(url, "Major aerospace companies (Google News)")
        candidates.sort(key=lambda item: item["date"], reverse=True)
        for aliases in organizations:
            match = next(
                (
                    item
                    for item in candidates
                    if item["url"] not in seen_urls
                    and not excluded_sources.search(item["source"])
                    and any(
                        re.search(r"\b" + re.escape(alias) + r"\b", item["title"] + " " + item["source"], re.IGNORECASE)
                        for alias in aliases
                    )
                ),
                None,
            )
            if match:
                selected.append(match)
                seen_urls.add(match["url"])
    selected.sort(key=lambda item: item["date"], reverse=True)
    return selected[:8]


def esa_items() -> list[dict[str, str]]:
    feeds = (
        ("https://www.esa.int/rssfeed/Our_Activities/Space_Science", "ESA Space Science"),
        ("https://www.esa.int/rssfeed/Our_Activities/Space_Transportation", "ESA Space Transportation"),
    )
    items = []
    for url, source in feeds:
        try:
            matching = rss_items(url, source)
            matching.sort(key=lambda item: item["date"], reverse=True)
            items.extend(matching[:2])
        except (URLError, TimeoutError, ET.ParseError, OSError) as error:
            print(f"Could not refresh {source}: {error}", file=sys.stderr)
    return deduplicate(items)


def major_space_news_items() -> list[dict[str, str]]:
    query = (
        '("space industry" OR "commercial space" OR satellite industry OR '
        'launch provider OR spaceflight OR space exploration) '
        '(launch OR satellite OR spacecraft OR mission OR constellation)'
    )
    regions = (("US", "US:en"), ("GB", "GB:en"), ("CA", "CA:en"), ("AU", "AU:en"), ("IN", "IN:en"), ("JP", "JP:en"))
    excluded_sources = re.compile(
        r"britannica|encyclopedia|motley fool|yahoo finance|investorplace|"
        r"stock titan|morningstar|marketwatch|business insider|tradingview|"
        r"investing\.com|finance|kalkine|dars\.gov\.et|guru.?focus|ad hoc news|stock|ticker",
        re.IGNORECASE,
    )
    relevant_title = re.compile(
        r"space|satellite|rocket|launch|spacecraft|orbital|constellation|"
        r"spaceflight|space station|lunar|mars",
        re.IGNORECASE,
    )
    items = []
    for region, edition in regions:
        url = "https://news.google.com/rss/search?" + urlencode(
            {"q": query, "hl": f"en-{region}", "gl": region, "ceid": edition}
        )
        try:
            items.extend(rss_items(url, "Major space news (Google News)"))
        except (URLError, TimeoutError, ET.ParseError, OSError) as error:
            print(f"Could not refresh major space news ({region}): {error}", file=sys.stderr)
    items = [item for item in items if not excluded_sources.search(item["source"]) and relevant_title.search(item["title"])]
    items.sort(key=lambda item: item["date"], reverse=True)
    return deduplicate(items)[:8]


def compute_tech_news_items() -> list[dict[str, str]]:
    queries = (
        (
            '(NVIDIA OR GPU OR supercomputer) ("fluid simulation" OR CFD OR turbulence)',
            "GPU-accelerated science and HPC",
        ),
        (
            'GPU accelerated computational fluid dynamics news',
            "GPU-accelerated science and HPC",
        ),
        (
            '("quantum computing" OR "quantum algorithm") ("fluid mechanics" OR "fluid dynamics" OR simulation)',
            "Quantum computing for fluid science",
        ),
        (
            '("quantum encoding" OR "quantum simulation") (fluids OR "fluid dynamics" OR CFD)',
            "Quantum computing for fluid science",
        ),
    )
    excluded_sources = re.compile(r"stock|finance|investing|motley fool|marketwatch|tradingview|guru.?focus|ad hoc news|pluang|wimi|award|billion|gaming|houdini|CG Channel|80 Level", re.I)
    technical_title = re.compile(r"GPU|NVIDIA|supercomput|fluid|CFD|turbulence|quantum|simulation|computational", re.I)
    gpu_title = re.compile(r"GPU|NVIDIA|supercomput|HPC|CFD|turbulence|computational fluid", re.I)
    items = []
    for query, topic in queries:
        url = "https://news.google.com/rss/search?" + urlencode(
            {"q": query, "hl": "en-US", "gl": "US", "ceid": "US:en"}
        )
        try:
            for item in rss_items(url, "GPU and quantum computing news"):
                item["topics"] = [topic]
                items.append(item)
        except (URLError, TimeoutError, ET.ParseError, OSError) as error:
            print(f"Could not refresh GPU/quantum news: {error}", file=sys.stderr)
    items = [
        item for item in items
        if not excluded_sources.search(item["source"] + " " + item["title"])
        and technical_title.search(item["title"])
        and (item["topics"][0] != "GPU-accelerated science and HPC" or gpu_title.search(item["title"]))
    ]
    items.sort(key=lambda item: item["date"], reverse=True)
    return deduplicate(items)[:8]


def deduplicate(items: list[dict[str, str]]) -> list[dict[str, str]]:
    seen: set[str] = set()
    unique = []
    for item in items:
        if item["url"] not in seen:
            seen.add(item["url"])
            unique.append(item)
    return unique


def industry_items() -> list[dict[str, str]]:
    aerospace_terms = re.compile(
        r"aerospace|aeronaut|aviation|aircraft|airline|rocket|launch|satellite|spacecraft|"
        r"commercial space|hypersonic|propulsion|air force|defense|defence|space industry|"
        r"space station|orbital|boeing|airbus|lockheed|spacex|blue origin",
        re.IGNORECASE,
    )
    by_source: dict[str, list[dict[str, str]]] = {}
    collectors = (
        ("SpaceX", spacex_items, INDUSTRY_ITEMS_PER_SOURCE),
        ("Blue Origin", blue_origin_items, INDUSTRY_ITEMS_PER_SOURCE),
        ("SpaceNews", lambda: rss_items("https://spacenews.com/feed/", "SpaceNews"), INDUSTRY_ITEMS_PER_SOURCE),
        ("NASA", lambda: rss_items("https://www.nasa.gov/news-release/feed/", "NASA"), INDUSTRY_ITEMS_PER_SOURCE),
        ("Major aerospace companies", major_company_items, 8),
        ("European Space Agency", esa_items, 4),
        ("Major space developments", major_space_news_items, 8),
        ("GPU and quantum computing", compute_tech_news_items, 6),
    )
    for source, collect, item_limit in collectors:
        try:
            matching = collect()
            if source in {"SpaceNews", "NASA"}:
                matching = [item for item in matching if aerospace_terms.search(item["title"] + " " + item["summary"][:500])]
            matching.sort(key=lambda item: item["date"], reverse=True)
            by_source[source] = matching[:item_limit]
        except (URLError, TimeoutError, ET.ParseError, OSError) as error:
            print(f"Could not refresh {source}: {error}", file=sys.stderr)
    items = [item for group in by_source.values() for item in group]
    items.sort(key=lambda item: item["date"], reverse=True)
    return deduplicate(items)[:sum(item_limit for _, _, item_limit in collectors)]


def spacewerx_sbirsttr_items() -> list[dict[str, str]]:
    html = fetch("https://spacewerx.us/get-funded/", {"User-Agent": "Mozilla/5.0"}).decode("utf-8", "replace")
    items = []
    call_pattern = re.compile(
        r'<a[^>]+href="([^"]+)"[^>]*>(DOW (?:SBIR|STTR) SPECIFIC TOPIC [^<]+)</a>(.{0,1400})',
        re.IGNORECASE | re.DOTALL,
    )
    today = datetime.now(timezone.utc).date()
    for match in call_pattern.finditer(html):
        title = plain_text(match.group(2))
        detail = plain_text(match.group(3))
        close_match = re.search(r"CLOSE\s*:?\s*(\d{1,2}\s+[A-Z]{3}\s+\d{2})", detail, re.I)
        open_match = re.search(r"OPEN\s*:?\s*(\d{1,2}\s+[A-Z]{3}\s+\d{2})", detail, re.I)
        if not close_match:
            continue
        try:
            close_date = datetime.strptime(close_match.group(1).title(), "%d %b %y").date()
            open_date = datetime.strptime(open_match.group(1).title(), "%d %b %y").date() if open_match else None
        except ValueError:
            continue
        if close_date < today or (open_date and open_date > today):
            continue
        items.append({
            "title": title,
            "url": "https://spacewerx.us/get-funded/",
            "date": close_date.strftime("%m/%d/%Y"),
            "source": "SpaceWERX / Department of War",
            "status": f"Open · closes {close_date.strftime('%b %d, %Y')}",
            "summary": "Current Space Force small-business solicitation. See the official SpaceWERX page for the topic description and submission portal.",
            "topics": ["Aerospace SBIR/STTR", "Space Force"],
        })
    return items


def sbir_gov_topic_items() -> list[dict[str, str]]:
    parser = SbirTopicListParser()
    try:
        parser.feed(fetch("https://www.sbir.gov/topics").decode("utf-8", "replace"))
    except (URLError, TimeoutError, OSError) as error:
        print(f"Could not refresh SBIR.gov active topics: {error}", file=sys.stderr)
        return []

    relevant = re.compile(
        r"aerospace|space|rocket|launch|satellite|aircraft|aviation|hypersonic|"
        r"propulsion|fluid|flow|turbulence|rarefied|multiphase|quantum|GPU|laser|"
        r"unmanned aerial|UAS",
        re.IGNORECASE,
    )
    results = []
    today = datetime.now(timezone.utc).date()
    for path, title in parser.topics:
        try:
            detail = plain_text(fetch("https://www.sbir.gov" + path).decode("utf-8", "replace"))
        except (URLError, TimeoutError, OSError) as error:
            print(f"Could not inspect SBIR.gov topic {path}: {error}", file=sys.stderr)
            continue
        status_match = re.search(r"Solicitation Status\s*:?\s*(Open|Closed|Archived)", detail, re.I)
        if not status_match or status_match.group(1).lower() != "open":
            continue

        description_match = re.search(
            r"\bDescription\s+(.*?)(?=\s+Site Map\s+Privacy Policy\b|$)",
            detail,
            re.I,
        )
        description = description_match.group(1).strip() if description_match else ""
        scope = f"{title} {description}"
        relevant_scope = re.compile(
            r"aerospace|aeronautic|spaceflight|spacecraft|satellite|rocket|launch|"
            r"hypersonic|high.speed flow|fluid|flow|turbulence|rarefied|propulsion|"
            r"aircraft|aviation|unmanned aerial|\bUAS\b",
            re.I,
        )
        if not relevant_scope.search(scope):
            continue

        date_format = r"(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)[a-z]*\s+\d{1,2},\s+\d{4}"
        open_match = re.search(rf"Open Date\s*:?\s*({date_format})", detail, re.I)
        close_match = re.search(rf"Close Date\s*:?\s*({date_format})", detail, re.I)
        closing = ""
        try:
            open_date = datetime.strptime(open_match.group(1), "%B %d, %Y").date() if open_match else None
            close_date = datetime.strptime(close_match.group(1), "%B %d, %Y").date() if close_match else None
        except ValueError:
            continue
        if open_date and open_date > today:
            continue
        if close_date:
            if close_date < today:
                continue
            closing = close_date.strftime("%m/%d/%Y")
        tags = [name for name, pattern in {
            "Space/aerospace": re.compile(r"aerospace|aeronautic|spaceflight|spacecraft|satellite|rocket|launch|aircraft|aviation|uas", re.I),
            "Hypersonics": re.compile(r"hypersonic|high.speed flow", re.I),
            "Fluids": re.compile(r"fluid|flow|turbulence|rarefied|multiphase|propulsion", re.I),
        }.items() if pattern.search(scope)]
        if re.search(r"Tagged as:.*?\b(?:SBIR|STTR|BOTH)\b", detail, re.I | re.S):
            tags.append("SBIR/STTR")
        results.append({
            "title": title,
            "url": "https://www.sbir.gov" + path,
            "date": closing,
            "source": "SBIR.gov active topics",
            "status": "Open · deadline not posted" if not closing else f"Open · closes {closing}",
            "summary": description[:420],
            "topics": tags,
        })
        time.sleep(0.1)
    return results


def nasa_sbirsttr_items() -> list[dict[str, str]]:
    try:
        html = fetch("https://www.nasa.gov/sbir_sttr/phase-ii-sequential/").decode("utf-8", "replace")
    except (URLError, TimeoutError, OSError) as error:
        print(f"Could not refresh NASA SBIR/STTR calls: {error}", file=sys.stderr)
        return []
    items = []
    anchor_pattern = re.compile(r'<a\b[^>]*href="(https://sam\.gov/workspace/contract/opp/[^"]+)"[^>]*>(.*?)</a>', re.I | re.S)
    for match in anchor_pattern.finditer(html):
        label = plain_text(match.group(2))
        if "2026 Phase II Sequential Call for White Papers" not in label:
            continue
        items.append({
            "title": "NASA 2026 SBIR Phase II Sequential Call for White Papers",
            "url": safe_url(match.group(1)),
            "date": "",
            "source": "NASA SBIR/STTR",
            "status": "Call listed · verify current submission window on SAM.gov",
            "summary": "NASA’s official SBIR/STTR page lists this Phase II Sequential call for eligible NASA Phase II awardees. The linked SAM.gov notice is authoritative for current dates and submission status.",
            "topics": ["Aerospace SBIR/STTR", "NASA"],
        })
    return deduplicate(items)


def funding_items() -> list[dict[str, str]]:
    candidates: dict[str, dict[str, str]] = {}
    title_relevance = re.compile(
        r"aerospace|spaceflight|space technology|space science|space grant|space exploration|"
        r"satellite|orbital|rocket|launch|lunar|mars|heliophysics|hypersonic|rarefied|non[- ]?equilibrium|kinetic theory|dsmc|"
        r"multiphase|multi[- ]phase|two[- ]phase|fluid|turbulence|aerodynamic|"
        r"quantum|chemistry|chemical|combustion|plasma|artificial intelligence|"
        r"machine learning|deep learning|computational fluid|transport phenomena|"
        r"propulsion|aircraft|aviation|rocket|launch vehicle|air force|space force|"
        r"sbir|sttr|small business innovation research|small business technology transfer|"
        r"gpu|high.performance computing|national laborator",
        re.IGNORECASE,
    )
    small_business_title = re.compile(
        r"\bSBIR\b|\bSTTR\b|small business innovation research|small business technology transfer",
        re.IGNORECASE,
    )
    small_business_agency = re.compile(
        r"NASA|National Aeronautics|Air Force|Space Force|Army|Navy|Naval Research|"
        r"Defense|\bDOE\b|Energy|National Science Foundation",
        re.IGNORECASE,
    )
    aerospace_sbirsttr_agency = re.compile(
        r"NASA|National Aeronautics|Air Force|Space Force|Army|Navy|Naval Research|Defense|\bDOD\b",
        re.IGNORECASE,
    )
    umbrella_title = re.compile(
        r"research interests|basic research|broad agency announcement|\bBAA\b|"
        r"transport phenomena|chemical process systems|space technology research",
        re.IGNORECASE,
    )
    search_terms = (
        "space", "spaceflight", "aerospace", "hypersonic", "rarefied flow", "rarefied gas",
        "multiphase flow", "two-phase flow", "quantum", "chemistry", "chemical kinetics",
        "artificial intelligence", "machine learning", "fluid dynamics", "flow physics",
        "turbulence", "computational fluid dynamics", "combustion", "plasma", "propulsion",
        "aerospace SBIR", "aerospace STTR", "NASA SBIR", "NASA STTR", "Air Force SBIR",
        "Space Force SBIR", "Army SBIR aerospace", "hypersonic SBIR", "propulsion SBIR",
        "Space Force research", "Army research laboratory aerospace", "national laboratory fluid dynamics",
        "Army STTR aerospace", "NASA small business aerospace", "national laboratory aerospace SBIR",
        "DOE national laboratory aerospace", "GPU scientific computing", "GPU research computing",
        "quantum computing fluid dynamics", "DOE laboratory fluid research",
    )
    endpoint = "https://api.grants.gov/v1/api/search2"
    for keyword in search_terms:
        start_record = 0
        hit_count = 1
        while start_record < hit_count:
            payload = json.dumps(
                {
                    "keyword": keyword,
                    "oppStatuses": "forecasted|posted",
                    "rows": 100,
                    "startRecordNum": start_record,
                }
            ).encode()
            request = Request(
                endpoint,
                data=payload,
                headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
            )
            try:
                with urlopen(request, timeout=25) as response:
                    data = json.load(response).get("data", {})
            except (URLError, TimeoutError, OSError, json.JSONDecodeError) as error:
                print(f"Could not refresh Grants.gov ({keyword}, {start_record}): {error}", file=sys.stderr)
                break

            hit_count = int(data.get("hitCount", 0))
            results = data.get("oppHits", [])
            if not results:
                break
            for opportunity in results:
                identifier = str(opportunity.get("id", ""))
                title = plain_text(str(opportunity.get("title", "")))
                status = str(opportunity.get("oppStatus", "")).lower()
                agency = plain_text(str(opportunity.get("agency", "")))
                is_sbirsttr = bool(small_business_title.search(title))
                if (
                    not identifier
                    or not title
                    or (is_sbirsttr and not small_business_agency.search(agency))
                    or not title_relevance.search(title)
                    and not umbrella_title.search(title)
                    or status not in {"posted", "forecasted"}
                ):
                    continue
                candidate = candidates.setdefault(identifier, dict(opportunity))
                if len(title) > len(str(candidate.get("title", ""))):
                    candidate.update(opportunity)
            start_record += len(results)
            time.sleep(0.15)

    topic_patterns = {
        "space": re.compile(r"aerospace|spaceflight|space technology|space science|space grant|space exploration|satellite|orbital|rocket|launch|lunar|mars|heliophysics|astrophysics", re.I),
        "hypersonic": re.compile(r"hypersonic|high[- ]speed flow", re.I),
        "rarefied": re.compile(r"rarefied|rarefaction|non[- ]equilibrium|kinetic theory|dsmc|knudsen", re.I),
        "multiphase": re.compile(r"multiphase|multi[- ]phase|two[- ]phase|particle[- ]laden|droplets?|bubbles?|cavitation", re.I),
        "quantum": re.compile(r"quantum", re.I),
        "chemistry": re.compile(r"chemistry|chemical kinetics|chemical process|combustion|reactive flows?", re.I),
        "flow physics": re.compile(r"fluid dynamics|fluid mechanics|flow physics|turbulence|aerodynamics|computational fluid|transport phenomena|navier.stokes|reacting flows?", re.I),
        "GPU/HPC computing": re.compile(r"GPU|graphics processing|high.performance computing|supercomputing|\bHPC\b", re.I),
    }
    ai_pattern = re.compile(r"artificial intelligence|machine learning|deep learning|physics-informed|\bAI\b", re.I)
    scientific_ai = re.compile(r"scientific discovery|science and energy|research infrastructure|science infrastructure", re.I)
    irrelevant_medical = re.compile(r"clinical trial|patient|human subjects|nervous system|mental health|cancer|opioid|diabetes|drug products|cybersecurity education|cyberai|scholarship for service", re.I)
    obsolete_title = re.compile(r"\b20\d{2}\s*[\u2012-\u2015-]\s*20\d{2}\b")
    today = datetime.now(timezone.utc).date()
    results = []
    for identifier, opportunity in candidates.items():
        title = plain_text(str(opportunity.get("title", "")))
        status = str(opportunity.get("oppStatus", "")).lower()
        closing = str(opportunity.get("closeDate", ""))
        year_range = obsolete_title.search(title)
        if year_range:
            years = re.findall(r"20\d{2}", year_range.group())
            if years and int(years[-1]) < today.year:
                continue
        try:
            if closing and datetime.strptime(closing, "%m/%d/%Y").date() < today:
                continue
        except ValueError:
            continue

        detail = {}
        needs_detail = not closing or bool(umbrella_title.search(title)) or bool(small_business_title.search(title))
        if needs_detail:
            try:
                request = Request(
                    "https://api.grants.gov/v1/api/fetchOpportunity",
                    data=json.dumps({"opportunityId": int(identifier)}).encode(),
                    headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
                )
                with urlopen(request, timeout=20) as response:
                    detail = json.load(response).get("data", {})
                time.sleep(0.1)
            except (URLError, TimeoutError, OSError, ValueError, json.JSONDecodeError) as error:
                print(f"Could not inspect Grants.gov opportunity {identifier}: {error}", file=sys.stderr)

        synopsis = detail.get("synopsis", {}) if isinstance(detail, dict) else {}
        detail_text = plain_text(str(synopsis.get("synopsisDesc", "")))
        scope = f"{title} {detail_text}"
        response_date = plain_text(str(synopsis.get("responseDateDesc", "")))
        archive_date = str(synopsis.get("archiveDate", ""))
        if re.search(r"archived|proposals are not accepted|no longer accepting", response_date, re.I):
            continue
        if archive_date:
            try:
                archive = datetime.strptime(archive_date[:12], "%b %d, %Y").date()
                if archive < today:
                    continue
            except ValueError:
                pass

        matches_topic = any(pattern.search(scope) for pattern in topic_patterns.values())
        has_ai_flow_scope = bool(ai_pattern.search(scope) and topic_patterns["flow physics"].search(scope))
        has_scientific_ai_scope = bool(ai_pattern.search(scope) and scientific_ai.search(scope))
        is_aerospace_sbirsttr = bool(small_business_title.search(title) and aerospace_sbirsttr_agency.search(str(opportunity.get("agency", ""))))
        if not matches_topic and not has_ai_flow_scope and not has_scientific_ai_scope and not is_aerospace_sbirsttr:
            continue
        if irrelevant_medical.search(title):
            continue

        if not closing and response_date:
            closing_match = re.search(r"\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)[a-z]*\.?\s+\d{1,2},?\s+\d{4}\b", response_date, re.I)
            if closing_match:
                try:
                    closing = datetime.strptime(closing_match.group().replace(".", ""), "%b %d, %Y").strftime("%m/%d/%Y")
                    if datetime.strptime(closing, "%m/%d/%Y").date() < today:
                        continue
                except ValueError:
                    closing = ""

        matched_topics = [name for name, pattern in topic_patterns.items() if pattern.search(scope)]
        if is_aerospace_sbirsttr:
            matched_topics.append("Aerospace SBIR/STTR")
        if has_ai_flow_scope:
            matched_topics.append("AI for flow physics")
        if has_scientific_ai_scope:
            matched_topics.append("AI for scientific research")
        sentences = re.split(r"(?<=[.!?])\s+", detail_text)
        relevant_sentences = [
            sentence for sentence in sentences
            if any(topic_patterns[name].search(sentence) for name in matched_topics if name in topic_patterns)
            or ("AI for flow physics" in matched_topics and ai_pattern.search(sentence))
            or ("AI for scientific research" in matched_topics and scientific_ai.search(sentence))
        ]
        summary = " ".join((relevant_sentences or sentences)[:2])[:420]
        results.append({
            "title": title,
            "url": f"https://www.grants.gov/search-results-detail/{quote(identifier)}",
            "date": closing,
            "source": plain_text(str(opportunity.get("agency", "Grants.gov"))),
            "status": (
                "Forecast · deadline not posted" if status == "forecasted" and not closing
                else "Forecast" if status == "forecasted"
                else "Open · deadline not posted" if not closing
                else "Open"
            ),
            "summary": summary,
            "topics": matched_topics,
        })
    results.sort(key=lambda item: (item["date"] or "9999", item["title"]))
    results.extend(spacewerx_sbirsttr_items())
    results.extend(sbir_gov_topic_items())
    results.extend(nasa_sbirsttr_items())
    results.extend(foundation_items())
    unique_results = []
    seen_opportunities: set[tuple[str, str]] = set()
    for item in results:
        key = (item["url"], item["title"])
        if key not in seen_opportunities:
            seen_opportunities.add(key)
            unique_results.append(item)
    results = unique_results
    results.sort(key=lambda item: (item.get("date") or "9999", item["title"]))
    return results


def foundation_items() -> list[dict[str, str]]:
    today = datetime.now(timezone.utc).date()
    opportunities = [
        {
            "title": "Simons Foundation Targeted Grants in Mathematics and the Physical Sciences",
            "url": "https://www.simonsfoundation.org/grant/targeted-grants-in-mps/",
            "date": "",
            "source": "Simons Foundation",
            "status": "Rolling · LOI required",
            "summary": (
                "Rolling letters of intent support high-risk theoretical mathematics, physics, "
                "and computer science. This is a broad physical-sciences program; confirm that "
                "a proposed fluid-mechanics project fits with the foundation before applying."
            ),
            "topics": ["Broad physical sciences; fluid fit to confirm"],
        },
        {
            "title": "Simons Collaborations in Mathematics and the Physical Sciences",
            "url": "https://www.simonsfoundation.org/grant/simons-collaborations-in-mathematics-and-the-physical-sciences/",
            "date": "10/29/2026",
            "source": "Simons Foundation",
            "status": "Open · LOI deadline",
            "summary": (
                "The foundation invites letters of intent for collaborative research in mathematics "
                "and the physical sciences. Program themes vary; confirm the current theme supports "
                "the proposed fluid-science topic."
            ),
            "topics": ["Broad physical sciences; fluid fit to confirm"],
        },
    ]
    active = []
    for opportunity in opportunities:
        if opportunity["date"] and datetime.strptime(opportunity["date"], "%m/%d/%Y").date() < today:
            continue
        active.append(opportunity)
    return active


def recent_industry_items(items: list[dict[str, str]]) -> list[dict[str, str]]:
    cutoff = datetime.now(timezone.utc).date() - timedelta(days=14)
    recent = []
    for item in items:
        try:
            if datetime.strptime(item["date"], "%Y-%m-%d").date() >= cutoff:
                recent.append(item)
        except (KeyError, TypeError, ValueError):
            continue
    return recent


def read_previous() -> dict:
    for path in (SOURCE_OUTPUT, SITE_OUTPUT):
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
    return {"industry": [], "papers": [], "funding": []}


def main() -> None:
    previous = read_previous()
    sections = {"industry": industry_items, "papers": paper_items, "funding": funding_items}
    data = {
        "updated": datetime.now(timezone.utc).isoformat(),
        "sources": ["SpaceX", "Blue Origin", "SpaceNews", "NASA", "ESA", "major aerospace companies", "GPU and quantum computing", "arXiv", "AIAA", "Journal of Fluid Mechanics", "Physics of Fluids", "Physical Review Fluids", "Grants.gov", "SpaceWERX", "SBIR.gov", "Simons Foundation"],
    }
    for section, collect in sections.items():
        try:
            fresh = collect()
            data[section] = fresh if fresh else previous.get(section, [])
        except (URLError, TimeoutError, ET.ParseError, OSError, json.JSONDecodeError) as error:
            print(f"Could not refresh {section}: {error}", file=sys.stderr)
            data[section] = previous.get(section, [])
        if section == "industry":
            data[section] = recent_industry_items(data[section])
    serialized = json.dumps(data, ensure_ascii=True, indent=2) + "\n"
    for path in (SOURCE_OUTPUT, SITE_OUTPUT):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(serialized, encoding="utf-8")
    print("Updated feed: " + ", ".join(f"{key}={len(data[key])}" for key in sections))


if __name__ == "__main__":
    main()