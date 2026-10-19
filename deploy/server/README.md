# Self-hosting on a Linux server (Docker + Caddy)

Runs the app on a server with a dedicated IP address and domain, with automatic HTTPS.
Works on x86_64 and ARM64 servers (e.g. 2 vCPU, 4 GB RAM, Ubuntu 24.04, plus 2–4 GB swap). The app needs about 3 GB RAM.

* Incoming: ports 80/443 are published only on the address in `PDAC_IP`; the Shiny container publishes no port.
* Outgoing: traffic from the containers (data download, certificate renewal) leaves with source address `PDAC_IP`.

```bash
# 1. Docker and Compose (Ubuntu archive packages)
sudo apt update && sudo apt install -y docker-compose-v2 git

# 2. Code and configuration
sudo git clone https://github.com/ANadern-2026/pdac-tissue-explorer.git /opt/pdac-tissue-explorer
cd /opt/pdac-tissue-explorer/deploy/server
sudo cp .env.example .env && sudo nano .env        # set PDAC_IP, PDAC_DOMAIN, ACME_EMAIL

# 3. Build, download the data (~245 MB from Zenodo, checksum-verified), start
sudo docker compose build
sudo docker compose run --rm app Rscript get_data.R
sudo docker compose up -d
```

Update to a new version: `cd /opt/pdac-tissue-explorer && sudo git pull && cd deploy/server && sudo docker compose up -d --build`.
Logs: `sudo docker compose logs -f app` (or `caddy`).
