from datetime import datetime, timezone
import os

from pymongo import MongoClient

MONGO_URI = os.environ.get("MONGODB_URI") or os.environ.get("MONGO_URI")
client = MongoClient(MONGO_URI)
db = client["nmf_tool"]
jobs_col = db["cnmf_jobs"]


def create_job(sc_job_id: str, email: str):
    jobs_col.insert_one({
        "sc_job_id": sc_job_id,
        "email": email,
        "status": "running",
        "created_at": datetime.now(timezone.utc),
        "completed_at": None,
        "s3_result_key": None,
        "error": None,
    })


def update_job(sc_job_id: str, status: str, s3_result_key: str = None, error: str = None):
    update = {
        "status": status,
        "completed_at": datetime.now(timezone.utc),
    }
    if s3_result_key:
        update["s3_result_key"] = s3_result_key
    if error:
        update["error"] = error
    jobs_col.update_one({"sc_job_id": sc_job_id}, {"$set": update})


def get_job(sc_job_id: str):
    return jobs_col.find_one({"sc_job_id": sc_job_id}, {"_id": 0})
