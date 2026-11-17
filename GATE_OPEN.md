# Gate G3: one line left on the droplet

G4 is closed. The repo is public at https://github.com/bellcheddar/KINFOLK, `main` is pushed
(d6f8c86) and the live snapshot is on `live` (8a7721c, 3,258 files, 164 MB).

## What is needed from Marc

In the DigitalOcean web console for the droplet, as root, one line:

    curl -fsSL https://raw.githubusercontent.com/bellcheddar/KINFOLK/main/deploy/pull/setup.sh | bash

It installs the service, the nginx site, the pull timer and the TLS certificate, and ends
`KINFOLK_READY`. DNS already resolves kinfolk.mdeller.com to the droplet.

## Then

Nothing by hand. The droplet's `kinfolk-pull.timer` fetches `live` every minute and publishes it,
so later releases need only `bash deploy/deploy.sh` on the Mac. Tell Claude when the console line
has finished and it will check every page with `deploy/check_live.sh`.

## State of the build

- Phases 0 to 10 complete. Phase 4 finished this session: the pocket trees and the Plate III
  tanglegram build (303 and 518 leaves), after neighbour joining was taught to drop domains
  whose KLIFS fingerprints have too few comparable rows.
- QC: 20 pipeline tests, 16 app tests, 16 of 16 page checks with no failures, axe clean,
  all five science checks pass.
- Outstanding: the one fold batch that failed on MASTL (821 residues) is folding again; when it
  lands the atlas is rebuilt and redeployed, which the pull timer picks up on its own.
- Pre-publication scan: no credentials. The droplet IP is in two earlier commits; it is public
  through DNS, so the history was left intact.
