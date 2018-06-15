/* Gaia DR2 analysis

https://www.ast.cam.ac.uk/ioa/wikis/WSDB/index.php/WSDB_GAIADR2

gaia_dr2.gaia_source
gaia_dr2_aux.gaia_source_neighbour_count

wsdb=> select * from gaia_dr2_aux.gaia_source_neighbour_count limit 1;
      source_id      |        ra        |        dec        | cnt0_25 | cnt0_50 | cnt1 | cnt2 | cnt4 | cnt8 | cnt16 | cnt32 | cnt64 | cnt128
 4164666032716788352 | 268.158775207367 | -8.79686399880187 |       1 |       1 |    1 |    1 |    1 |    4 |    10 |    36 |   151 |    584


explore the contents

MOD(random_index, 10000) = 12 

*/


/* Gaia triplets */
/* could include outer join with WISE and PS */
/* do some stats too */
/* mult cnt128 x 252 to get sources density per square degree */

/*
*/

/* b>89; 80 sec;     7943 rows */
/* b>88; 82 sec;    32269 rows */
/* b>85; 84 sec;   204373 rows */
/* b>80; 90 sec;   817775 rows */
/* b>75; 102 sec; 1925624 rows */
/* b>60; 748 sec; 8063081 rows */
/* b>45;  sec; 204373 rows */
/* 29.5<b<30.5;  124 secs;  3086613 rows */
/* 29.5<b<30.5 
   AND 90 < l < 270 ;  40sec;  970835 rows 
   AND 90 > l > 270 ;  104sec; 2115778 rows */
/* 24.5<b<25.5;  sec;  rows */
/* 19.75>b>20.25;  139 sec; 3552408 rows */
/* b>15;  sec; 204373 rows */
/* b>10;  sec; 204373 rows */
/* b>5;  sec; 204373 rows */
/* b>0;  sec; 204373 rows */
/* all; */ 

SELECT
    CAST(0 AS INT) AS cnt4, count(*),
    MIN(cnt4), MAX(cnt4), MEDIAN(cnt4),
    MIN(cnt8), MAX(cnt8), MEDIAN(cnt8),
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile99,
    MIN(cnt64), MAX(cnt64), MEDIAN(cnt64),    
    MIN(cnt128), MAX(cnt128), MEDIAN(cnt128),
    MIN(cnt128 - cnt64), MAX(cnt128 - cnt64), MEDIAN(cnt128 - cnt64),
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile99

FROM
    /* gaia_dr2_aux.gaia_source_neighbour_count AS gsnc */
    gaia_dr2.gaia_source AS gs
    INNER JOIN gaia_dr2_aux.gaia_source_neighbour_count AS gsnc
        ON gs.source_id = gsnc.source_id 
WHERE
    /* b > 80 is around 1% of sky area; cf sin(80 degrees) */
    /* b > 80.0 */
    /* b > 48.951 is 25% of sky or 5,000 */
    /* b > 50.0 */
    /* 20 deg x 20 deg = 1% of sky */
    /* (gsnc.ra between 170.0 and 190.0)
    AND (gsnc.dec between -10.0 and 10.0) */
    b > 89.0
    AND l > 0.0
    AND l < 360.0 
    /* AND (b BETWEEN 29.75 AND 30.25)*/
    /* AND (l < 90.0 OR l > 270.0) */
    /* (b between 19.75 and 20.25) */
    /* (b between 14.95 and 15.05) */
    /* b > 60.0 */    
    /* (b < 15.0 and b > -15.0) */
    AND cnt4 >= 3

/*

UNION

SELECT
    CAST(1 AS INT) AS cnt4, count(*),
    MIN(cnt4), MAX(cnt4), MEDIAN(cnt4),
    MIN(cnt8), MAX(cnt8), MEDIAN(cnt8),
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt8) AS cnt8_entile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile99,
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile99

FROM
    gaia_dr2.gaia_source AS gs
    INNER JOIN gaia_dr2_aux.gaia_source_neighbour_count AS gsnc
        ON gs.source_id = gsnc.source_id
WHERE
    /* b > 80 is around 1% of sky area; cf sin(80 degrees) */
    /* b > 48.951 is 25% of sky or 5,000 */
    /* b > 50.0 */
    b > 60.0
    /* b > 60.0 */    
    /* (b < 15.0 and b > -15.0) */
    AND cnt4 = 1

UNION

SELECT
    CAST(2 AS INT) AS cnt4, count(*),
    MIN(cnt4), MAX(cnt4), MEDIAN(cnt4),    
    MIN(cnt8), MAX(cnt8), MEDIAN(cnt8),
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile99,
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile99

FROM
    gaia_dr2.gaia_source AS gs
    INNER JOIN gaia_dr2_aux.gaia_source_neighbour_count AS gsnc
        ON gs.source_id = gsnc.source_id
WHERE
    /* b > 80 is around 1% of sky area; cf sin(80 degrees) */
    /* b > 48.951 is 25% of sky or 5,000 */
    /* b > 50.0 */
    b > 60.0
    /* b > 60.0 */    
    /* (b < 15.0 and b > -15.0) */
    AND cnt4 = 2

UNION

SELECT
    CAST(3 AS INT) AS cnt4, count(*),
    MIN(cnt4), MAX(cnt4), MEDIAN(cnt4),
    MIN(cnt8), MAX(cnt8), MEDIAN(cnt8),
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile99,
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile99


FROM
    gaia_dr2.gaia_source AS gs
    INNER JOIN gaia_dr2_aux.gaia_source_neighbour_count AS gsnc
        ON gs.source_id = gsnc.source_id
WHERE
    /* b > 80 is around 1% of sky area; cf sin(80 degrees) */
    /* b > 48.951 is 25% of sky or 5,000 */
    /* b > 50.0 */
    b > 60.0
    /* b > 60.0 */    
    /* (b < 15.0 and b > -15.0) */
    AND cnt4 = 3

UNION

SELECT
    CAST(4 AS INT) AS cnt4, count(*),
    MIN(cnt4), MAX(cnt4), MEDIAN(cnt4),
    MIN(cnt8), MAX(cnt8), MEDIAN(cnt8),
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile99,
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile99

FROM
    gaia_dr2.gaia_source AS gs
    INNER JOIN gaia_dr2_aux.gaia_source_neighbour_count AS gsnc
        ON gs.source_id = gsnc.source_id
WHERE
    /* b > 80 is around 1% of sky area; cf sin(80 degrees) */
    /* b > 48.951 is 25% of sky or 5,000 */
    /* b > 50.0 */
    b > 60.0
    /* b > 60.0 */    
    /* (b < 15.0 and b > -15.0) */
    AND cnt4 = 4

UNION

SELECT
    CAST(5 AS INT) AS cnt4, count(*),
    MIN(cnt4), MAX(cnt4), MEDIAN(cnt4),
    MIN(cnt8), MAX(cnt8), MEDIAN(cnt8),
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile99,
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile99

FROM
    gaia_dr2.gaia_source AS gs
    INNER JOIN gaia_dr2_aux.gaia_source_neighbour_count AS gsnc
        ON gs.source_id = gsnc.source_id
WHERE
    /* b > 80 is around 1% of sky area; cf sin(80 degrees) */
    /* b > 48.951 is 25% of sky or 5,000 */
    /* b > 50.0 */
    b > 60.0
    /* b > 60.0 */    
    /* (b < 15.0 and b > -15.0) */
    AND cnt4 = 5


UNION

SELECT
    CAST(6 AS INT) AS cnt4, count(*),
    MIN(cnt4), MAX(cnt4), MEDIAN(cnt4),
    MIN(cnt8), MAX(cnt8), MEDIAN(cnt8),
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt8) AS cnt8_centile99,
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt128 - cnt64)
        AS cnt128m64_centile99

FROM
    gaia_dr2.gaia_source AS gs
    INNER JOIN gaia_dr2_aux.gaia_source_neighbour_count AS gsnc
        ON gs.source_id = gsnc.source_id
WHERE
    /* b > 80 is around 1% of sky area; cf sin(80 degrees) */
    /* b > 48.951 is 25% of sky or 5,000 */
    /* b > 50.0 */
    b > 60.0
    /* b > 60.0 */    
    /* (b < 15.0 and b > -15.0) */
    AND cnt4 > 5

ORDER BY cnt4 ASC

*/

/*
SELECT
    count(*),
    MIN(cnt8), MAX(cnt8), MEDIAN(cnt8),
    PERCENTILE_DISC(array[0.01,0.05,0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99])
        WITHIN GROUP (ORDER BY cnt8) AS centiles
FROM
    gaia_dr2.gaia_source AS gs
    INNER JOIN gaia_dr2_aux.gaia_source_neighbour_count AS gsnc
        ON gs.source_id = gsnc.source_id
WHERE
    /* b > 80 is around 1% of sky area; cf sin(80 degrees) */
    /* b > 48.951 is 25% of sky or 5,000 */
    /* b > 50.0 */
    b > 80.0
    /* b > 60.0 */ 
    /* (b < 15.0 and b > -15.0) */
    AND cnt4 > 2


UNION


SELECT
    CAST(1 AS INT) AS cnt4, count(*),
    MIN(cnt4), MAX(cnt4), MEDIAN(cnt4),
    percentile_disc(array[0.01,0.05,0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99])) WITHIN GROUP (ORDER BY cnt4) AS centiles,
FROM
    gaia_dr2.gaia_source AS gs
    INNER JOIN gaia_dr2_aux.gaia_source_neighbour_count AS gsnc
        ON gs.source_id = gsnc.source_id
WHERE
    /* b > 80 is around 1% of sky area; cf sin(80 degrees) */
    /* b > 48.951 is 25% of sky or 5,000 */
    /* b > 50.0 */
    b > 80.0
    /* b > 60.0 */    
    /* (b < 15.0 and b > -15.0) */
    cnt4 = 1


UNION

SELECT
    CAST(2 AS INT) AS cnt4, count(*),
    MIN(cnt4), MAX(cnt4), MEDIAN(cnt4),
    percentile_disc(0.01) WITHIN GROUP (ORDER BY cnt4) AS centile1,
    percentile_disc(0.05) WITHIN GROUP (ORDER BY cnt4) AS centile5,
    percentile_disc(0.10) WITHIN GROUP (ORDER BY cnt4) AS centile10,
    percentile_disc(0.25) WITHIN GROUP (ORDER BY cnt4) AS centile25,
    percentile_disc(0.50) WITHIN GROUP (ORDER BY cnt4) AS centile50,
    percentile_disc(0.75) WITHIN GROUP (ORDER BY cnt4) AS centile75,
    percentile_disc(0.90) WITHIN GROUP (ORDER BY cnt4) AS centile90,
    percentile_disc(0.95) WITHIN GROUP (ORDER BY cnt4) AS centile95,
    percentile_disc(0.99) WITHIN GROUP (ORDER BY cnt4) AS centile99
FROM
    gaia_dr2.gaia_source AS gs
    INNER JOIN gaia_dr2_aux.gaia_source_neighbour_count AS gsnc
        ON gs.source_id = gsnc.source_id
WHERE
    /* b > 80 is around 1% of sky area; cf sin(80 degrees) */
    /* b > 48.951 is 25% of sky or 5,000 */
    /* b > 50.0 */
    b > 80.0
    /* b > 60.0 */    
    /* (b < 15.0 and b > -15.0) */
    cnt4 = 2



ORDER BY cnt4

*/
