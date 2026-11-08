#!/usr/bin/env python3
"""
validate_profile.py — валидатор контракта профиля (profile v2).

Проверяет пару profile.tsv + profile.json (sidecar) на соответствие
docs/profile_contract.md: структура, типы, инварианты
(input = classified + unclassified; status <=> число строк;
согласованность идентичности TSV <-> sidecar; сумма reads <= classified).

Запускается:
  * в пайплайне — framework-процессом VALIDATE_PROFILE на каждый выход
    адаптера (нарушение контракта = падение задачи);
  * в юнит-тестах и ручной проверке адаптеров.

Опционально сверяет идентичность с ожидаемой (--expect-*).

Выход: 0 при успехе (печатает OK), 1 с сообщением при нарушении.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from profile_io import ProfileError, read_profile  # noqa: E402


def main():
    ap = argparse.ArgumentParser(
        description="Валидация пары profile.tsv + profile.json (profile v2)")
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--expect-sample", default=None)
    ap.add_argument("--expect-fraction", default=None)
    ap.add_argument("--expect-role", default=None)
    args = ap.parse_args()

    try:
        p = read_profile(args.tsv, args.json)
    except ProfileError as e:
        sys.exit(f"INVALID: {e}")

    for opt, key in (("--expect-sample", "sample_id"),
                     ("--expect-fraction", "fraction"),
                     ("--expect-role", "role")):
        want = getattr(args, opt[2:].replace("-", "_"))
        if want is not None and p[key] != want:
            sys.exit(f"INVALID: {args.tsv}: {key}={p[key]!r}, "
                     f"ожидалось {want!r} ({opt})")

    print(f"OK: {args.tsv} + {args.json} "
          f"({p['sample_id']}/{p['fraction']}/{p['role']}, "
          f"status={p['status']}, taxa={len(p['taxa'])}, "
          f"input_reads={p['input_reads']})")


if __name__ == "__main__":
    main()
