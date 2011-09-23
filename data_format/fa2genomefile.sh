 seqstat -a $1 | grep '^*' | awk '{print $2,$3}'
