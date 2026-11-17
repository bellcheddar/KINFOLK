# No gate is open

G3 and G4 both closed on 2026-11-17. https://kinfolk.mdeller.com is live and the repo is
public at https://github.com/bellcheddar/KINFOLK.

## How a release goes out now

    bash deploy/deploy.sh

It runs the tests, pushes the commit to `main` and a snapshot of `app/` plus the data bundle to
`live`, then waits for the droplet to publish it and checks every page. The droplet's
`kinfolk-pull.timer` fetches `live` every minute; nothing inbound over SSH is needed. The last
release went live about 50 seconds after the push.

One-time droplet setup, should it ever be needed again (it is idempotent):

    curl -fsSL https://raw.githubusercontent.com/bellcheddar/KINFOLK/main/deploy/pull/setup.sh | bash

Run it straight after a push and pin it, because raw.githubusercontent.com serves a branch path
from cache for a few minutes and will otherwise install the previous files:

    curl -fsSL https://raw.githubusercontent.com/bellcheddar/KINFOLK/<sha>/deploy/pull/setup.sh \
      | KINFOLK_REF=<sha> bash

## What was settled during deployment

- The droplet already runs other apps and one of them holds `127.0.0.1:8010`, so KINFOLK binds
  `/run/kinfolk/kinfolk.sock` instead. Until that was found, nginx served chatMCD at
  kinfolk.mdeller.com while `/static` came off disk and looked correct.
- TLS is in place for kinfolk.mdeller.com and renews through certbot's scheduled task.

## State of the build

- Phases 0 to 10 complete. Phase 4 finished in this session: the pocket trees and the Plate III
  tanglegram build (303 and 518 leaves). Plate VII's Q5 figure was added, so every question on
  the Almanac now carries one.
- QC: 20 pipeline tests, 16 app tests, 16 of 16 page checks with no failures, axe clean, all five
  science checks pass, locally and against the droplet.
- Outstanding: the fold batch that failed on MASTL (821 residues) is folding again on this Mac.
  When it lands the atlas is rebuilt and redeployed; nothing on the droplet needs touching.
