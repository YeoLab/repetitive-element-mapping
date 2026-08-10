# Beads Advanced: Full bd CLI Reference

`bd prime` covers basic CRUD. This file covers the advanced features you need for forge pipeline integration and sophisticated workflow management.

## Scriptable Output

- `bd --json <any-command>` returns machine-readable JSON on ANY command
- `bd q "Title" --type task` quick-capture: creates issue, outputs ONLY the ID (for scripting)
- `bd -q` (quiet mode) suppresses non-essential output

## Labels

```bash
bd label add <id> <label>           # Add label (e.g., stage:3-implement, blocker, review-finding)
bd label remove <id> <label>        # Remove label
bd label list <id>                  # List labels on an issue
bd label list-all                   # All unique labels in the database
bd label propagate <parent> <label> # Propagate label to all children
bd list -l <label> --flat           # Filter issues by label
```

Label conventions for forge: `phase-0` through `phase-7`, `stage:<stage-id>`, `pipeline-run`, `blocker`, `review-finding`, `security-finding`, `hitl`, `discovered-in:<stage>`.

## Issue Types

Beyond task/bug/feature, bd supports: `epic`, `decision` (ADR), `spike` (timeboxed investigation), `story`, `milestone`, `chore`. Use `--type=decision` for architectural decisions, `--type=spike` for research.

## Gates (Async Coordination)

Gates block workflow until a condition is met. Critical for HITL (human-in-the-loop) approval.

```bash
bd gate list                        # Show all open gates
bd gate list --all                  # Include closed gates
bd gate resolve <id>                # Close a gate (human approval)
bd gate check                       # Evaluate all open gates
bd gate show <id>                   # Gate details
bd human <id>                       # Flag issue for human decision
bd human list                       # Show all issues awaiting human input
bd human respond <id> "response"    # Respond to a human gate
bd human dismiss <id>               # Dismiss a human flag
```

Gate types: `human` (manual resolve), `timer` (auto-expire), `gh:run` (GitHub Actions), `gh:pr` (PR merge), `bead` (cross-rig bead closure).

## Formulas and Molecules (Workflow Templates)

Formulas are YAML workflow templates. Molecules are running instances.

```bash
# Formulas (templates)
bd formula list                     # Available templates
bd formula show <name>              # Show formula details

# Molecules (instances)
bd mol pour <formula> [--var key=val]   # Create persistent molecule from formula
bd mol wisp <formula> [--var key=val]   # Create ephemeral molecule (ops, releases)
bd mol progress <id>                # Show molecule completion progress
bd mol current                      # Show current position in active workflow
bd mol show <id>                    # Show molecule structure
bd mol bond <a> <b>                 # Combine two molecules/protos
bd mol bond <a> <b> --type parallel # Run in parallel
bd mol stale                        # Find complete-but-unclosed molecules
bd mol ready                        # Molecules ready for gate-resume dispatch
```

Formula search paths: `.beads/formulas/` (project), `~/.beads/formulas/` (user).

## State Dimensions (Operational Tracking)

Track operational state on issues using dimension:value labels, with event history.

```bash
bd set-state <id> <dim>=<val> --reason "why"   # Set state atomically (creates event + updates label)
bd state <id> <dim>                             # Query current value of a dimension
bd state list <id>                              # List all state dimensions
```

Example: `bd set-state forge-pipeline stage=3-implement --reason "Stage starting"` tracks pipeline position.

## Key-Value Store

Persistent key-value pairs for flags, counters, config.

```bash
bd kv set <key> <value>             # Store value
bd kv get <key>                     # Retrieve value
bd kv list                          # List all keys
bd kv clear <key>                   # Delete key
```

Example: `bd kv set cost:3-implement "2.45"` stores per-stage cost data.

## Query Language

Complex filtering beyond simple list flags.

```bash
bd query "status=open AND priority<=2"
bd query "(status=open OR status=blocked) AND priority<2"
bd query "type=bug AND label=urgent"
bd query "assignee=none AND type=task"
bd query "status=open AND updated>7d"         # Updated in last 7 days
```

## Batch Operations

```bash
bd create --graph plan.json         # Create graph of issues with deps from JSON
bd create --file tasks.md           # Batch create from markdown
bd close <id1> <id2> <id3> ...      # Close multiple at once
bd update <id1> <id2> --claim       # Claim multiple at once
```

Graph JSON format:
```json
{
  "issues": [
    {"id": "a", "title": "Task A", "type": "task"},
    {"id": "b", "title": "Task B", "type": "task", "deps": ["a"]}
  ]
}
```

## Dependency Graph Visualization

```bash
bd graph <id>                       # Terminal DAG
bd graph --compact <id>             # Tree format, scannable
bd graph --dot <id> | dot -Tsvg > g.svg   # SVG export
bd dep tree <id>                    # Dependency tree
bd dep add <issue> <depends-on>     # Add dependency (issue blocked by depends-on)
bd dep remove <issue> <depends-on>  # Remove dependency
```

## Epics and Hierarchy

```bash
bd epic create "Epic title" --description="..."   # Create epic
bd children <parent-id>             # List child beads of a parent
bd create --title="Subtask" --parent=<epic-id>     # Create child issue
```

## Swarms (Multi-Agent Coordination)

For parallel work across multiple agents on an epic.

```bash
bd swarm validate <epic-id>         # Check deps, find cycles, estimate parallelism
bd swarm create <epic-id>           # Create swarm for parallel execution
bd swarm status                     # Monitor swarm progress
```

## Merge-Slots (Serialized Conflict Resolution)

Prevents multiple agents from racing to resolve the same conflicts.

```bash
bd merge-slot create                # Create for current rig
bd merge-slot check                 # Is slot available?
bd merge-slot acquire               # Try to acquire
bd merge-slot release               # Release
```

## Export/Import

```bash
bd export -o backup.jsonl           # Export all issues + memories
bd export --no-memories -o issues.jsonl  # Issues only
bd import backup.jsonl              # Import (upsert semantics)
bd import --dry-run backup.jsonl    # Preview only
```

## Diagnostics

```bash
bd doctor                           # Full health check
bd doctor --agent --json            # Agent-friendly diagnostics with remediation commands
bd preflight                        # PR readiness checklist
bd stale                            # Issues with no recent activity
bd orphans                          # Issues with broken dependencies
bd lint                             # Check issues for missing template sections
```

## Valid Statuses

`open`, `in_progress`, `blocked`, `deferred`, `closed`, `pinned`, `hooked`

## Memories (Cross-Session Knowledge)

```bash
bd remember "insight text"          # Store persistent memory
bd remember "key insight" --key mykey  # Store with specific key
bd memories                         # List all memories
bd memories <keyword>               # Search memories
bd recall <key>                     # Retrieve specific memory
bd forget <key>                     # Remove memory
```

## Worktrees (Parallel Development)

```bash
bd worktree create <name>           # Create git worktree with shared beads db
bd worktree list                    # List active worktrees
bd worktree remove <name>           # Clean up worktree
```
