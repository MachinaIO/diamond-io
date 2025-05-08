#!/bin/sh

set -eux

EXIT_CODE=0
RUNNER_TAG=$1

SUBNET=$(./check_availability.sh ${INSTANCE_TYPE})

cat cloud-init.sh | sed -e "s#__GH_REPO__#${GH_REPO}#" -e "s/__GH_PAT__/${GH_PAT}/" -e "s/__RUNNER_TAG__/${RUNNER_TAG}/" > .startup.sh

INSTANCE_ID=$(aws ec2 run-instances \
  --image-id ${IMAGE_ID} \
  --block-device-mapping "[ { \"DeviceName\": \"/dev/sda1\", \"Ebs\": { \"VolumeSize\": 100, \"DeleteOnTermination\": true } } ]" \
  --ebs-optimized \
  --instance-initiated-shutdown-behavior terminate \
  --instance-type ${INSTANCE_TYPE} \
  --key-name devopsoregon \
  --security-group-ids ${SG} \
  --subnet-id ${SUBNET} \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=machina-io-ephemeral-${RUNNER_TAG}},{Key=ProjectName,Value=machina-io}]" "ResourceType=volume,Tags=[{Key=ProjectName,Value=machina-io}]" \
  --user-data "file://.startup.sh" \
  --query "Instances[0].InstanceId" \
  --output text) || { EXIT_CODE=1; exit $EXIT_CODE; }
 echo "INSTANCE_ID=$INSTANCE_ID" >> $GITHUB_ENV
 aws ec2 wait instance-running --instance-ids $INSTANCE_ID
 echo "EC2 instance $INSTANCE_ID is running"

exit $EXIT_CODE
