#!/bin/bash

# cancelruns.sh -- cancel runs on SLURM without all the faff.
# John Wilkinson 18/12/2020

ALL=0

function cancel {
  # ask to cancel run, and if yes, cancel it
  # pass this function a line in squeue -u $USER.
  # $1 = JOBID
  # $2 = PARTITION
  # $3 = NAME
  # $4 = USER
  # ... but not all these are strictly necessary!

  # if this is the title line of squeue, skip over it
  if [ "$1" == 'JOBID' ]; then
    return 0
  fi

  # if the user has previously asked to cancel all, bypass the prompt
  if [ $ALL == '1' ]; then
    scancel $1
    return 0
  fi

  # ask the user if they want to cancel (the < /dev/tty is important!)
  read -p "Cancel run $1 ($3)? (y(es)/n(o)/a(ll)/m(ore)/q(uit))" CHOICE < /dev/tty

  case "$CHOICE" in
    a|A ) echo "Cancelling all runs"; ALL=1; scancel $1;;
    y|Y ) scancel $1;;
    n|N ) echo "Not cancelling $1";;
    m|M ) squeue | grep $1; cancel $1;;
    q|Q ) exit;;
    * ) echo "Invalid choice, asking agian..."; cancel $1;
  esac
  return 0
}

squeue -u $USER | while read line; do cancel $line; done
