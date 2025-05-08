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

  IS_AVAILABLE=$(aws ec2 describe-instance-type-offerings --region "$REGION" --location-type availability-zone --filters Name=location,Values="$AZ" Name=instance-type,Values="$INSTANCE_TYPE" --query "InstanceTypeOfferings[].InstanceType" --output text 2>/dev/null)
  [ -n "$IS_AVAILABLE" ] && { SELECTED_SUBNET="$SUBNET_ID"; break; }
done

if [ -n "$SELECTED_SUBNET" ]; then
  echo "$SELECTED_SUBNET"
else
  echo "Error: No subnet supports $INSTANCE_TYPE"
  exit 1
fi
