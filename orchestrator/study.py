
import argparse
import optuna
from evaluate import evaluate
from schemas import Metrics

SPACE = {
  "cfl": [0.3, 0.5, 0.8],
  "rk_stages": [2,3],
  "compiler_flags": ["-O2","-O3","-O3 -ffast-math"]
}

def objective(trial):
    cfg = {k: trial.suggest_categorical(k, v) for k,v in SPACE.items()}
    m: Metrics = evaluate(cfg)
    trial.set_user_attr("metrics", m.model_dump())
    return m.score

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--trials", type=int, default=5)
    args = p.parse_args()
    study = optuna.create_study(direction="maximize")
    study.optimize(objective, n_trials=args.trials)
    print("Best score:", study.best_trial.value)
    print("Best params:", study.best_trial.params)
    print("Best metrics:", study.best_trial.user_attrs.get("metrics"))
