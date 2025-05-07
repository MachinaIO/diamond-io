#!/bin/sh

set -eux
trap 'poweroff' TERM EXIT INT

GH_REPO="__GH_REPO__"
GH_PAT="__GH_PAT__"
RUNNER_TAG="__RUNNER_TAG__"

# idle poweroff script
cat >> /etc/crontab << 'EOF'
*/1 * * * * root cat /proc/uptime | awk -F ' ' '{ if ($1 < 300) exit 1 }' && (cat /proc/loadavg | awk -F ' ' '{ if ($1 <= .3 && $2 < .3 && $3 < .3) exit 1 }' || poweroff)
EOF
GH_RUNNER_TOKEN=$(curl -s -X POST -H "Authorization: Bearer ${GH_PAT}" -H "Accept: application/vnd.github.v3+json" https://api.github.com/repos/${GH_REPO}/actions/runners/registration-token | jq -r '.token')

# basics
apt-get update
apt-get install -y curl jq git ca-certificates gnupg lsb-release

# github runner
mkdir actions-runner && cd actions-runner
curl -o actions-runner-linux-x64-2.323.0.tar.gz -L https://github.com/actions/runner/releases/download/v2.323.0/actions-runner-linux-x64-2.323.0.tar.gz
tar xzf ./actions-runner-linux-x64-2.323.0.tar.gz

cd /
chown -R ubuntu actions-runner
cd actions-runner
sudo -u ubuntu ./config.sh --url https://github.com/${GH_REPO} --token ${GH_RUNNER_TOKEN} --labels ${RUNNER_TAG} --ephemeral --unattended --disableupdate --replace
sudo -u ubuntu ./run.sh
