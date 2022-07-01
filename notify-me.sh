#!/bin/bash
datdir=$1
program=$2
if [ "$#" -ne 2 ]; then
echo " Usage: cmd [datdir] [program]"
exit 1
fi

set -eu

# Slack Webhook IncomingのAPIへの通知メッセージ用のJSONを作る。
gen_post_data()
{
  cat <<EOF
{
    "blocks": [
        {
            "type": "section",
            "text": {
                "type": "mrkdwn",
                "text": "dataset ${datdir} finish \n program ${program}"
            }
         }
    ]
}
EOF
}

# 作ったメッセージを付けてSlack Webhook IncomingのAPIを叩く。
curl -s -S -X POST -d "$(gen_post_data)" "https://hooks.slack.com/services/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"