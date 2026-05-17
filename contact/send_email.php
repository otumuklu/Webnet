<?php
if ($_SERVER["REQUEST_METHOD"] == "POST") {
    $name = trim($_POST['name'] ?? '');
    $email = trim($_POST['email'] ?? '');
    $subject = trim($_POST['subject'] ?? '');
    $message = trim($_POST['message'] ?? '');

    if ($name === '' || $subject === '' || $message === '' || !filter_var($email, FILTER_VALIDATE_EMAIL)) {
        header("Location: index.html?status=error");
        exit();
    }

    $to = "tumuko@rpi.edu";
    $email_subject = "HAVA contact form: " . str_replace(["\r", "\n"], '', $subject);
    $email_message = "Name: " . $name . "\n";
    $email_message .= "Email: " . $email . "\n\n";
    $email_message .= "Message:\n" . $message . "\n";

    $headers = "From: HAVA Website <tumuko@rpi.edu>\r\n" .
               "Reply-To: " . $email . "\r\n" .
               "X-Mailer: PHP/" . phpversion();

    if (mail($to, $email_subject, $email_message, $headers)) {
        header("Location: index.html?status=success");
    } else {
        header("Location: index.html?status=error");
    }
    exit();
}

header("Location: index.html");
exit();
?>
