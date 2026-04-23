import boto3

ses = boto3.client("ses", region_name="us-east-2")
SENDER = "azahanbabaev@gmail.com"


def send_job_complete(recipient: str, sc_job_id: str):
    ses.send_email(
        Source=SENDER,
        Destination={"ToAddresses": [recipient]},
        Message={
            "Subject": {"Data": "Your cNMF Analysis is Complete"},
            "Body": {
                "Text": {
                    "Data": (
                        f"Your single-cell cNMF job has finished.\n\n"
                        f"Job ID: {sc_job_id}\n\n"
                        f"Return to the NMF Tool and enter your Job ID in the "
                        f"'Retrieve Previous Job' section to load your results."
                        f"You have 24 hours to access your results before they are deleted."
                    )
                }
            },
        },
    )


def send_job_failed(recipient: str, sc_job_id: str):
    ses.send_email(
        Source=SENDER,
        Destination={"ToAddresses": [recipient]},
        Message={
            "Subject": {"Data": "Your cNMF Analysis Failed"},
            "Body": {
                "Text": {
                    "Data": (
                        f"Unfortunately your single-cell cNMF job failed.\n\n"
                        f"Job ID: {sc_job_id}\n\n"
                        f"Please return to the NMF Tool to retry or contact support."
                        f"You can see the error by accessing previous job results on the tool page."
                    )
                }
            },
        },
    )
