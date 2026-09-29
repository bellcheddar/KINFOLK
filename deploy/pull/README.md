# Pull deploys (HTTPS only from the Mac)

`bash deploy/deploy.sh` pushes the code to `main` and a snapshot of `app/` plus the data bundle
to `live`. The droplet's `kinfolk-pull.timer` publishes `live` within a minute and writes
`/static/version.txt`. No SSH from the Mac; the repo is public, so the droplet needs no key.

## One-time setup

In the DigitalOcean web console for the droplet, as root, type:

```
curl -fsSL https://raw.githubusercontent.com/bellcheddar/KINFOLK/main/deploy/pull/setup.sh | bash
```

It creates the `kinfolk` service user, installs `uv`, the systemd service, the nginx site, the
pull script and its timer, publishes `live` once, obtains the TLS certificate and ends with
`KINFOLK_READY`. It is idempotent: rerun it after any failure.

## Operating it

| Want | Do |
|---|---|
| Deploy | `bash deploy/deploy.sh` |
| See what is live | `curl https://kinfolk.mdeller.com/static/version.txt` |
| Check every page | `bash deploy/check_live.sh` |
| Droplet log | `/var/log/kinfolk-pull.log`, `journalctl -u kinfolk` |
| Stop auto-publishing | `systemctl disable --now kinfolk-pull.timer` |
| Roll back | `git push --force origin <old-live-commit>:refs/heads/live` |
| Emergency (SSH) | `bash deploy/deploy-rsync.sh`, then `bash deploy/repair_perms.sh` if pages return 500 |
