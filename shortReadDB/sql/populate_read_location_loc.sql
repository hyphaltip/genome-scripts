use sReadDB;
INSERT INTO read_location_loc (read_location_id, location) SELECT read_location_id, 
	LineFromText(CONCAT('LINESTRING(',chromosome_id, ' ',startmatch,',',chromosome_id, ' ',
	endmatch,')')) FROM read_location;
