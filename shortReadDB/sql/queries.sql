
-- count the features that overlap 5' UTRs
use sReadDB;

--EXPLAIN 
SELECT c.name, r.read_id, ri.machid, r.startmatch,r.endmatch, ri.trim_3,  
       f.feature_id, f.fid,f.fstart, f.fstop,f.fstrand
FROM read_location AS r, read_info ri, feature AS f, chromosome AS c
WHERE f.ftype = 'five_prime_utr' AND
      ri.trim_3 >= 18 AND ri.trim_3 <= 28 AND      
      c.chromosome_id = r.chromosome_id AND
      ri.read_id = r.read_id AND
      r.chromosome_id = f.chromosome_id AND
      MBRIntersects(f.location,r.location)
INTO OUTFILE '/var/tmp/five_prime_utr.dat';

SELECT c.name, r.read_id,ri.machid,r.startmatch,r.endmatch, ri.trim_3,  
       f.feature_id, f.fid, f.fstart, f.fstop,f.fstrand
FROM read_location AS r, read_info ri, feature AS f, chromosome AS c
WHERE f.ftype = 'three_prime_utr' AND
      ri.trim_3 >= 18 AND ri.trim_3 <= 28 AND      
      c.chromosome_id = r.chromosome_id AND
      ri.read_id = r.read_id AND
      r.chromosome_id = f.chromosome_id AND
      MBRIntersects(f.location,r.location)
INTO OUTFILE '/var/tmp/three_prime_utr.dat';

SELECT c.name, r.read_id, ri.machid,r.startmatch,r.endmatch, ri.trim_3,  
       f.feature_id, f.fid,f.fstart, f.fstop,f.fstrand
FROM read_location AS r, read_info ri, feature AS f, chromosome AS c
WHERE f.ftype = 'mrna' AND
      ri.trim_3 >= 18 AND ri.trim_3 <= 28 AND      
      c.chromosome_id = r.chromosome_id AND
      ri.read_id = r.read_id AND
      r.chromosome_id = f.chromosome_id AND
      MBRIntersects(f.location,r.location)
INTO OUTFILE '/var/tmp/mrna.dat';


/*
EXPLAIN
SELECT r.read_id, r.startmatch,r.endmatch,  
       f.feature_id, f.fstart, f.fstop
FROM read_location AS r, feature AS f
WHERE f.ftype = 'five_prime_utr' AND
      r.chromosome_id = f.chromosome_id AND
      r.chromosome_id = 1 AND
      MBRIntersects(r.location,f.location)
--      ORDER by r.chromosome_id, f.fstart
*/
