# Self-hosting on a Linux server (Docker + Caddy)

Runs the app on a server with a dedicated IP address and domain, with automatic HTTPS.
Works on x86_64 and ARM64 servers (e.g. 2 vCPU, 4 GB RAM, Ubuntu 24.04, plus 2–4 GB swap). The app reads each gene's values on demand from the on-disk data store and needs about 0.5 GB RAM; building the store once (step 3) needs about 2 GB.

* Incoming: ports 80/443 are published only on the address in `PDAC_IP`; the Shiny container publishes no port.
* Outgoing: traffic from the containers (data download, certificate renewal) leaves with source address `PDAC_IP`.

```bash
# 1. Docker and Compose (Ubuntu archive packages)
sudo apt update && sudo apt install -y docker-compose-v2 git

# 2. Code and configuration
sudo git clone https://github.com/ANadern-2026/pdac-tissue-explorer.git /opt/pdac-tissue-explorer
cd /opt/pdac-tissue-explorer/deploy/server
sudo cp .env.example .env && sudo nano .env        # set PDAC_IP, PDAC_DOMAIN, ACME_EMAIL

# 3. Build, download the data (~245 MB from Zenodo, checksum-verified) and build the data store (~10 min), start
sudo docker compose build
sudo docker compose run --rm app Rscript get_data.R
sudo docker compose up -d
```

Update to a new version: `cd /opt/pdac-tissue-explorer && sudo git pull && cd deploy/server && sudo docker compose up -d --build`.
If the release notes say the data changed, run `sudo docker compose run --rm app Rscript get_data.R` before `up -d` (it only downloads or rebuilds what is missing or outdated).
Logs: `sudo docker compose logs -f app` (or `caddy`).

## Usage counts

The app counts visits and in-app usage anonymously (`app/R/05_usage.R`): a random visit id, start and end time, the referring site's host name, a screen-size class, an optional `?ref=` tag in the address, and the pages, tabs, genes and downloads used. No cookies, IP addresses or browser details are stored. The counts are in `usage.sqlite` in the Docker volume `pdac-usage`, separate from the data volume.

* Set `PDAC_STATS_KEY` and `PDAC_STATS_PASSWORD` in `.env` (e.g. `openssl rand -hex 16` each) and run `sudo docker compose up -d`; then open `https://<domain>/?stats=<PDAC_STATS_KEY>` and enter the password. Without the exact key the app shows the normal site, with no prompt or hint of the page. The visit that opens the page is not counted, and the page can stop counting your own browser.
* Tag links you post with `?ref=...` (e.g. `https://<domain>/?ref=linkedin`) to see them separately.
* Back up the counts: `sudo docker compose cp app:/srv/pdac/usage/usage.sqlite ./usage-$(date +%F).sqlite`.
* Local runs (RStudio) count nothing: counting needs the variable `PDAC_USAGE_DB`, which only `compose.yaml` sets.
