<?php
if ($_SERVER["REQUEST_METHOD"] == "POST") {
    // Get form data
    $name = $_POST['firstname'] . ' ' . $_POST['lastname'];
    $email = $_POST['email'];
    $subject = $_POST['subject'];

    // Set recipient email address
    $to = 'tumuko@rpi.edu';

    // Set email subject
    $email_subject = 'New Contact Form Submission';

    // Build the email content
    $email_message = "Name: $name\n";
    $email_message .= "Email: $email\n\n";
    $email_message .= "Subject:\n$subject\n";

    // Set email headers
    $headers = "From: $name <$email>";

    // Attempt to send the email
    if (mail($to, $email_subject, $email_message, $headers)) {
        // If the email is sent successfully, redirect back to the contact form with a success message
        header("Location: contact.html?status=success");
        exit();
    } else {
        // If there is an error sending the email, redirect back to the contact form with an error message
        header("Location: contact.html?status=error");
        exit();
    }
} else {
    // If it's not a POST request, redirect back to the contact form
    header("Location: contact.html");
    exit();
}
?>
