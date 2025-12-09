
from pydantic import BaseModel

class Metrics(BaseModel):
    l2: float
    runtime_ms: float
    determinism: float
    score: float
