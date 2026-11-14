"""Turn published IceCube alerts into the neutrino events the lane consumes.

    python3 scripts/ingest_icecube.py                    # the bundled sample
    python3 scripts/ingest_icecube.py --input my.json    # your own archive
    python3 scripts/ingest_icecube.py --check            # validate, write nothing

Accepts the GCN Kafka schema (``ra_dec_error`` at 90% containment, ``nu_energy``
in TeV, ``far`` in Hz, ``trigger_time``) and the legacy GCN/AMON table layout
(``error90_arcmin``, ``energy_tev``, ``far_per_year``). Every conversion applied
to every row is recorded, because the difference between a 90% radius and a 1
sigma is a factor of 2.1 in the cone and a factor of ~4.5 in the chance
probability, which is larger than most of the effects this repository measures.

Rows are rejected, not repaired, and each rejection is written into
``results/icecube_ingest.json`` with its reason. A retracted alert, a row with
no trigger time, or an error with no stated containment is a reason to stop, not
to fill in a default.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from seti_search.icecube import IceCubeSchemaError, from_alert  # noqa: E402

DEFAULT_INPUT = os.path.join(ROOT, "data", "icecube_alerts_sample.json")
DEFAULT_OUT = os.path.join(ROOT, "data", "neutrino_events.json")
REPORT = os.path.join(ROOT, "results", "icecube_ingest.json")


def read_rows(source: str) -> tuple[list[dict], str]:
    """Read a JSON archive from a path or an http(s) URL."""
    if source.startswith(("http://", "https://")):
        try:
            with urllib.request.urlopen(source, timeout=20) as response:
                raw = response.read().decode("utf-8")
        except Exception as exc:  # noqa: BLE001 - reported, never hidden
            raise SystemExit(f"could not fetch {source}: {exc}")
    else:
        with open(source, encoding="utf-8") as fh:
            raw = fh.read()
    payload = json.loads(raw)
    if isinstance(payload, dict):
        for key in ("alerts", "events", "notices"):
            if key in payload:
                payload = payload[key]
                break
        else:
            payload = [payload]
    if not isinstance(payload, list):
        raise SystemExit("the archive is neither a list of alerts nor an object "
                         "containing one")
    return [row for row in payload if isinstance(row, dict)], source


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument("--allow-test", action="store_true",
                        help="accept hardware injections, labelled as test events")
    parser.add_argument("--check", action="store_true",
                        help="validate and report without writing the catalog")
    args = parser.parse_args()

    rows, source = read_rows(args.input)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    accepted, rejected = [], []
    for index, row in enumerate(rows):
        provenance = str(row.get("_source") or source)
        release = f"published IceCube alert ingested from {provenance}"
        try:
            parsed = from_alert(row, public_source_url=provenance,
                                data_release=release,
                                allow_test=args.allow_test)
        except (IceCubeSchemaError, ValueError, KeyError) as exc:
            rejected.append({
                "index": index,
                "event_name": row.get("event_name", row.get("id")),
                "reason": str(exc),
            })
            continue
        accepted.append(parsed)

    events = [a.event for a in accepted]
    ledger = {
        "source": source,
        "ingested_utc": stamp,
        "n_rows": len(rows),
        "n_accepted": len(accepted),
        "n_rejected": len(rejected),
        "rejected": rejected,
        "alerts": [a.to_dict() for a in accepted],
        "conversions": sorted({c for a in accepted for c in a.conversions}),
        "note": ("every accepted row is a published alert; the only quantities "
                 "changed are units and the containment-to-sigma conversion, "
                 "both recorded above"),
    }

    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    with open(REPORT, "w", encoding="utf-8") as fh:
        json.dump(ledger, fh, indent=2)

    if not args.check:
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump([e.to_dict() for e in events], fh, indent=2)

    print(f"ICECUBE INGEST  ({source})")
    print(f"   {'event':24} {'time UTC':22} {'RA':>8} {'Dec':>7} "
          f"{'sigma':>6} {'TeV':>9} {'grade':8} {'FAR/yr':>8}")
    for alert in accepted:
        e = alert.event
        far = "-" if alert.far_per_year is None else f"{alert.far_per_year:.3g}"
        print(f"   {e.event_id:24} {e.event_time_utc.strftime('%Y-%m-%dT%H:%M:%S'):22} "
              f"{e.ra_deg:8.2f} {e.dec_deg:7.2f} {e.angular_error_deg:6.3f} "
              f"{e.energy_estimate_gev / 1000.0:9.2f} {e.alert_status:8} {far:>8}")
    for row in rejected:
        print(f"   REJECTED row {row['index']}: {row['reason'][:100]}")
    if accepted:
        print("\n   conversions applied:")
        for conversion in ledger["conversions"]:
            print(f"     {conversion}")
    print(f"\n   {len(accepted)} accepted, {len(rejected)} rejected")
    if args.check:
        print("   --check: no catalog written")
    else:
        print(f"   wrote {os.path.relpath(args.out, ROOT)} "
              "(this is now the lane's real catalog)")
    print(f"   wrote {os.path.relpath(REPORT, ROOT)}")
    return 1 if rejected else 0


if __name__ == "__main__":
    raise SystemExit(main())
