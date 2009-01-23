/* project */

-- DROP TABLE IF EXISTS project;

CREATE TABLE IF NOT EXISTS project (
	project_id	MEDIUMINT UNSIGNED NOT NULL AUTO_INCREMENT,
	name		varchar(32) NOT NULL,
	code		varchar(32) NULL, /* like 200M6AAXX_4 from BCGSC for Ncrass*/
	description	varchar(255) NULL, /* */
 	PRIMARY KEY (project_id),
	UNIQUE KEY ui_name (name),
	UNIQUE KEY ui_code (code),
	KEY i_desc (description)
);

/* short-read info */

--DROP TABLE IF EXISTS read_info;

CREATE TABLE IF NOT EXISTS read_info (
	read_id	INTEGER(10) UNSIGNED NOT NULL AUTO_INCREMENT,
	project_id MEDIUMINT UNSIGNED NOT NULL REFERENCES project(project_id),
	machid	varchar(32) NOT NULL, /* machine assigned Id for read */
	seq	char(36) NOT NULL,
	qual	char(36) NOT NULL DEFAULT '',
	adaptor_trimmed 	boolean NOT NULL DEFAULT FALSE,
	trim_5  		TINYINT NOT NULL DEFAULT '-1', /* 5' trim site */
	trim_3  		TINYINT NOT NULL DEFAULT '-1', /* 3' trim site */
	PRIMARY KEY (read_id),
	KEY        i_project (project_id),
	UNIQUE KEY ui_machid (project_id,machid),
	KEY        i_seq  (seq),
	KEY        i_qual (qual),
	KEY i_trimmed (adaptor_trimmed),
	KEY i_trim (trim_3,trim_5)
);

/* chromosome info */
-- DROP TABLE IF EXISTS chromosome;

CREATE TABLE IF NOT EXISTS chromosome (
	chromosome_id	MEDIUMINT UNSIGNED NOT NULL AUTO_INCREMENT,
	organism	varchar(255) NOT NULL,
	name		varchar(64) NOT NULL,
	length		MEDIUMINT UNSIGNED NOT NULL,
	description	varchar(255) NULL,
	PRIMARY KEY (chromosome_id),
	UNIQUE KEY ui_name_organism (organism,name),
	KEY i_length (length)
);


/* read mapped location info */
-- DROP TABLE IF EXISTS analysis_method;

CREATE TABLE IF NOT EXISTS analysis_method (
	method_id	SMALLINT UNSIGNED NOT NULL AUTO_INCREMENT,
	program		varchar(64) NOT NULL,
	program_params	varchar(255) NOT NULL DEFAULT '',
	
	PRIMARY KEY (method_id),
	UNIQUE KEY ui_prog (program,program_params)
);

/* read mapped location info */
-- DROP TABLE IF EXISTS read_location;
CREATE TABLE IF NOT EXISTS read_location (
	read_location_id	integer(10) UNSIGNED NOT NULL AUTO_INCREMENT,
	method_id		SMALLINT UNSIGNED NOT NULL REFERENCES analysis_method(method_id),
	chromosome_id		MEDIUMINT UNSIGNED NOT NULL REFERENCES chromosome(chromosome_id),
	read_id                 integer(10) UNSIGNED NOT NULL REFERENCES read_info (read_id),
	startmatch		MEDIUMINT UNSIGNED NOT NULL,
	endmatch		MEDIUMINT UNSIGNED NOT NULL,
	strand			boolean NOT NULL DEFAULT '1', /* 1 = +, 0 = - */
	score			TINYINT UNSIGNED NOT NULL DEFAULT '0',
	mismatches		TINYINT UNSIGNED NOT NULL DEFAULT '0',
	quality			TINYINT UNSIGNED NOT NULL DEFAULT '0', /* sum quality scores in trimmed seq? */	
	PRIMARY KEY (read_location_id),
	KEY i_read (read_id),
	KEY i_method (method_id),
	KEY i_chrom (chromosome_id),
	KEY i_loc (startmatch,endmatch),
	KEY i_end (endmatch),
	KEY i_strand (strand),
	KEY i_score (score),
	KEY i_mm (mismatches),
	KEY i_qual (quality)
);

-- DROP TABLE IF EXISTS read_location_loc;

CREATE TABLE IF NOT EXISTS read_location_loc (
	read_location_id	integer(10) UNSIGNED NOT NULL,
	location		LINESTRING NOT NULL,
	PRIMARY KEY (read_location_id),
	SPATIAL INDEX (location)

);
/* features */
--DROP TABLE IF EXISTS feature;

CREATE TABLE IF NOT EXISTS feature (
	feature_id		MEDIUMINT UNSIGNED NOT NULL AUTO_INCREMENT,
	chromosome_id		MEDIUMINT UNSIGNED NOT NULL REFERENCES chromosome(chromosome_id),
	fid			varchar(128) NOT NULL,
	fname			varchar(128) NULL,
	fparent			varchar(128) NULL,
	ftype			varchar(64) NOT NULL,
	fstart			MEDIUMINT UNSIGNED NOT NULL,
	fstop			MEDIUMINT UNSIGNED NOT NULL,
	fstrand			boolean NOT NULL DEFAULT '1', /* 1 = +, 0 = - */
	fsource			varchar(64) NOT NULL,
	location		LINESTRING NOT NULL,	
	PRIMARY KEY (feature_id),
        KEY i_chrom (chromosome_id),
	KEY ui_fid (fid),
	KEY i_name (fname),
        KEY i_type (ftype),
	KEY i_parent (fparent),
	KEY i_source (fsource),
	KEY i_start (fstart),
	KEY i_stop   (fstop),
	KEY i_strand (fstrand),
	SPATIAL INDEX (location)
);
