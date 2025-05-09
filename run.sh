#!/bin/sh

set -eux

EXIT_CODE=0
RUNNER_TAG=$1

VPC_ID="vpc-085ffb1026b00654e"
SUBNET_IDS=$(aws ec2 describe-subnets --filters "Name=vpc-id,Values=${VPC_ID}" "Name=tag:Type,Values=private" --query "Subnets[].SubnetId" --output json | jq -r '.[]')


cat cloud-init.sh | sed -e "s#__GH_REPO__#${GH_REPO}#" -e "s/__GH_PAT__/${GH_PAT}/" -e "s/__RUNNER_TAG__/${RUNNER_TAG}/" > .startup.sh

for SUBNET in $SUBNET_IDS; do
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
    --output text 2>&1 || true)
    if echo "$INSTANCE_ID" | grep -q '^i-[0-9a-f]\{17\}$'; then
      echo "INSTANCE_ID=$INSTANCE_ID" >> "$GITHUB_ENV"
      aws ec2 wait instance-running --instance-ids "$INSTANCE_ID"
      echo "EC2 instance $INSTANCE_ID is running in subnet $SUBNET"
      exit 0  # EXIT_CODE=0
    elif echo "$INSTANCE_ID" | grep -q "InsufficientInstanceCapacity"; then
      echo "Warning: Insufficient capacity for $INSTANCE_TYPE in subnet $SUBNET"
      EXIT_CODE=1
      continue
    else
      echo "Error: Failed to launch EC2 instance with subnet $SUBNET: $INSTANCE_ID"
      exit 1
    fi
done

echo -e "Error: No subnets in \n$SUBNET_IDS \nhad sufficient capacity for $INSTANCE_TYPE."
exit $EXIT_CODE
