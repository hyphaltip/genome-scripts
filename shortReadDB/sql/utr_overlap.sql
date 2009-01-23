select f.fid,b.fid from feature as f, feature as b WHERE f.ftype = 'three_prime_utr' AND b.ftype = 'three_prime_utr' AND f.chromosome_id = b.chromosome_id AND MBRIntersects(f.location,b.location) AND f.feature_id != b.feature_id AND f.fstrand != b.fstrand
INTO OUTFILE '/var/tmp/three_prime_utr.overlap_set.dat';
