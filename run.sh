#!/bin/sh

set -eux

EXIT_CODE=0
RUNNER_TAG=$1

cat cloud-init.sh | sed -e "s#__GH_REPO__#${GH_REPO}#" -e "s/__GH_PAT__/${GH_PAT}/" -e "s/__RUNNER_TAG__/${RUNNER_TAG}/" > .startup.sh

INSTANCE_ID=$(aws ec2 run-instances \
  --image-id ${IMAGE_ID} \
  --block-device-mapping "[ { \"DeviceName\": \"/dev/sda1\", \"Ebs\": { \"VolumeSize\": 16, \"DeleteOnTermination\": true } } ]" \
  --ebs-optimized \
  --instance-initiated-shutdown-behavior terminate \
  --instance-type ${INSTANCE_TYPE} \
  --key-name devops \
  --security-group-ids sg-0e511f05c162bb458 \
  --subnet-id subnet-0e1c657227636974b \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=drill-${RUNNER_TAG}},{Key=ProjectName,Value=devops}]" "ResourceType=volume,Tags=[{Key=ProjectName,Value=devops}]" \
  --user-data "file://.startup.sh" \
  --query "Instances[0].InstanceId" \
  --output text) || EXIT_CODE=1
 echo "INSTANCE_ID=$INSTANCE_ID" >> $GITHUB_ENV
 aws ec2 wait instance-running --instance-ids $INSTANCE_ID
 echo "EC2 instance $INSTANCE_ID is running"
