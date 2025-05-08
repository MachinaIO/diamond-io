#!/bin/bash
#
# Script to check if an EC2 instance type is available in a subnet's availability zone

REGION="us-west-2"
INSTANCE_TYPE="$1"
VPC_ID="vpc-085ffb1026b00654e"
SUBNET_IDS=$(aws ec2 describe-subnets --filters "Name=vpc-id,Values=${VPC_ID}" "Name=tag:Type,Values=private" --query "Subnets[].SubnetId" --output json | jq -r '.[]')
SELECTED_SUBNET=""

[ -z "$SUBNET_IDS" ] && { echo "Error: No private subnets found in VPC $VPC_ID"; exit 1; }

for SUBNET_ID in $SUBNET_IDS; do
  AZ=$(aws ec2 describe-subnets --subnet-ids "$SUBNET_ID" --region "$REGION" --query "Subnets[0].AvailabilityZone" --output text 2>/dev/null)
  [ -z "$AZ" ] && continue

  RESULT=$(aws ec2 run-instances \
    --region "$REGION" \
    --instance-type "$INSTANCE_TYPE" \
    --image-id "$IMAGE_ID" \
    --subnet-id "$SUBNET_ID" \
    --security-group-ids $SG \
    --tag-specifications "ResourceType=instance,Tags=[{Key=ProjectName,Value=machina-io}]" "ResourceType=volume,Tags=[{Key=ProjectName,Value=machina-io}]" \
    --count 1 \
    --dry-run \
    --query "Instances[0].InstanceId" \
    --output text 2>&1 || echo "ERROR: $?")

  if [[ "$RESULT" == *"Request would have succeeded, but DryRun flag is set"* ]]; then
    SELECTED_SUBNET="$SUBNET_ID"
    SELECTED_AZ="$AZ"
    break
  elif [[ "$RESULT" == *"InsufficientInstanceCapacity"* ]]; then
    echo "Insufficient capacity in $AZ"
  elif [[ "$RESULT" == *"InvalidParameterCombination"* ]]; then
    echo "Failed to check capacity in $AZ (Invalid parameter combination: $RESULT)"
  else
    echo "Failed to check capacity in $AZ (Error: $RESULT)"
  fi
done


if [ -n "$SELECTED_SUBNET" ]; then
  echo "$SELECTED_SUBNET"
else
  echo "Error: No subnet supports $INSTANCE_TYPE"
  exit 1
fi
