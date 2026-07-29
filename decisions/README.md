# Decision Log

This directory tracks every methodological decision made during the
project. Each decision is recorded as a standalone Markdown file in
the [Architecture Decision Record (ADR)](https://adr.github.io/)
format, adapted for a research project.

## Why we keep this log

Research is full of choices that look obvious in hindsight but were
non-trivial at the time. Three months later, a reviewer (or an
advisor) will ask "why did you pick X and not Y?" — and without a
record, the answer becomes guesswork.

Keeping this log:
- forces us to articulate *why* a choice was made, not just *that* it
  was made;
- protects the project against drift — we can revisit and revise a
  decision explicitly instead of silently changing direction;
- gives a future supervisor a single place to audit the project's
  reasoning.

## File naming

`NNNN-short-slug.md` where `NNNN` is a four-digit sequence number
starting from `0001`. The slug is a hyphen-separated phrase summarising
the decision.

## Format of each ADR

Each file has six sections:

1. **Status** — `proposed`, `accepted`, `superseded by NNNN`, or
   `deprecated`.
2. **Context** — what problem are we facing? What constraints apply?
3. **Decision** — what we chose.
4. **Alternatives considered** — what else was on the table, and why
   it was rejected.
5. **Consequences** — what becomes easier, what becomes harder, what
   we commit to as a result.
6. **Date and author** — when, and by whom.

## Style

Plain prose. No marketing language. State trade-offs honestly. If a
decision was made under uncertainty, say so; do not pretend it was
obvious.

## Index

| ID   | Title                                        | Status              |
|------|----------------------------------------------|---------------------|
| 0001 | Dataset choice for Phase 2                   | superseded by 0003  |
| 0002 | Quantum-walk conventions                     | accepted            |
| 0003 | Switch to short-term dataset                 | accepted            |
| 0004 | Seizure duration policy                      | accepted            |
| 0005 | Patient-level analysis                       | accepted            |
| 0006 | ID1 as development patient                   | accepted            |

New entries get appended to this table as they are added.