""""


TODO: look at the W1-W2 colours for all Gaia pairs; also look at WISE
nion-detections

https://www.ast.cam.ac.uk/ioa/wikis/WSDB/index.php/WSDB_GAIADR2
gaia_dr2.gaia_source

See here
https://www.ast.cam.ac.uk/ioa/wikis/WSDB/index.php/WSDB_GAIADR2_AUX#Neighbor_table

wsdb=> select * from gaia_dr2_aux.gaia_source_neighbour_count limit 5;

      source_id      |        ra        |        dec        | cnt0_25 | cnt0_50 | cnt1 | cnt2 | cnt4 | cnt8 | cnt16 | cnt32 | cnt64 | cnt128
---------------------+------------------+-------------------+---------+---------+------+------+------+------+-------+-------+-------+--------
 4164666032716788352 | 268.158775207367 | -8.79686399880187 |       1 |       1 |    1 |    1 |    1 |    4 |    10 |    36 |   151 |    584
 4164666032721252864 | 268.162921175373 |  -8.7992920387406 |       1 |       1 |    1 |    1 |    4 |    5 |    12 |    41 |   147 |    599
 4164666032716765056 |  268.16327175346 | -8.79868137404748 |       1 |       1 |    1 |    2 |    4 |    4 |    11 |    36 |   151 |    600
 4164666028420860160 | 268.163489504561 | -8.79893229296021 |       1 |       1 |    1 |    2 |    4 |    4 |    11 |    37 |   150 |    597
 4164666032721250816 | 268.163979711806 | -8.79930172774628 |       1 |       1 |    1 |    1 |    4 |    4 |    11 |    38 |   149 |    595
(5 rows)


Requires:
https://github.com/segasai/astrolibpy/tree/master/utils





Requires:
https://github.com/segasai/astrolibpy/tree/master/utils

this is the easyaccess git commit that adds the SQL to the FITS format
output file:

https://github.com/mgckind/easyaccess/commit/aa1ea0adfc059993ff5dfb72fd156237d5088b4d

fits[1].write_history('Created by easyaccess ' + version.__version__ + ' on ' + created)
fits[1].write_comment('Query = ' + query)


WISE columns are described here:

http://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec2_1a.html

https://www.ast.cam.ac.uk/ioa/wikis/WSDB/index.php/WSDB_WISE
ra            | double precision  |
dec           | double precision  |
w1mpro        | real              |
w1sigmpro     | real              |
w1snr         | real              |
w2mpro        | real              |
w2sigmpro     | real              |
w2snr         | real              |
w3mpro        | real              |
w3sigmpro     | real              |
w3snr         | real              |
w4mpro        | real              |
w4sigmpro     | real              |
w4snr         | real              |


# Aperture magnitudes could be useful for data screening
w1mag_1       | real              |
w1sigm_1      | real              |
w1flg_1      | integer           |
w2mag_1       | real              |
w2sigm_1      | real              |
w2flg_1       | integer           |
w3mag_1       | real              |
w3sigm_1      | real              |
w3flg_1       | integer           |
w4mag_1       | real              |
w4sigm_1      | real              |
w4flg_1       | integer           |

w1mcor W1 aperture curve-of-growth correction, in magnitudes. This correction
is subtracted from the nominal 8.25" aperture photometry brightness, w1mag_2,
to give the "standard-aperture" magnitude.

These look worth looking at too and a few aperture for W1 and W2; ignore
W3 and W4 aperture magniutdes for now




"""

from __future__ import print_function

# standard library functions
import time
t0 = time.time()

import os
import sys

import datetime
import argparse

import sqlparse

# 3rd party functions
# import atpy
# print(atpy.__version__)
from astropy.table import Table
from astropy.io import fits
import pandas as pd

# private functions
sys.path.append('/home/rgm/soft/python/lib/')
from librgm import sqlutil
# print(sqlutil.__version__)
# help(sqlutil)
from librgm.wsdb import rd_config_wsdb
from librgm.table_metadata import table_metadata

# print time to load all the python functions
print('Elapsed time(secs): ', time.time() - t0)



def table_add_sqlquery(table=None, query=None):
    """

    common indentation is removed with textwrap.dedent()

    import textwrap
    query = textwrap.dedent(query)

    """

    def cutquery(query, length):
        """
        Return query in a list of fixed sized character strings

        suggest length = 70 so can be stored in FITS format COMMENT
        from easyaccess codebase
        """

        result = [query[0 + i:length + i] for i in range(0, len(query), length)]

        return result

    table.meta['SQL'] = query
    comment = []
    comment.append('START:SQLQUERY')
    for iline, line in enumerate(query.splitlines()):
        print(iline, ':', line)
        cutresult = cutquery(line, 70)
        print('cutquery:', len(cutresult))
        for icut, record in enumerate(cutresult):
            comment.append(record)

    comment.append('END:SQLQUERY')
    table.meta['COMMENT'] = comment

    return table




def get_storedquery(queryname=None):
    """

    various stored queries

    could store as a dict:

    storedquery = {
            'query1': query1,
            'query2': query2,
    }


    storedquery['query3'] = query3

    or

    storedquery.update({'query3': query3})
    which can update with multiple key, value pairs



    common indentation is removed with textwrap.dedent()

    Should not SELECT on RA, Dec for both tables

    This is what Sergey says we should do:

    with x as (select * wise where ra>.. ra< w1< ... )
    y as  (select * from x, ps1 where q3c_join() )
    select * from z where rmag-w1<33


    An Example from Sergey's Wiki


    Get the SDSS and 2MASS photometry for the BHB candidates
    (according to Yanny 2000 criteria):

    SELECT s.ra, s.dec, s.psfmag_u, s.psfmag_g,
           s.psfmag_r, s.psfmag_i, s.psfmag_z, p.j_m, p.h_m, p.k_m
    FROM sdssdr7.phototag as s, twomass.psc AS p
    WHERE (psfmag_u-psfmag_g)<1.35 AND (psfmag_u-psfmag_g)>0.8
                AND (psfmag_g-psfmag_r)<0 AND (psfmag_g-psfmag_r)>-0.4
                AND type=6 and mode=1
                AND  q3c_join(s.ra, s.dec, p.ra, p.decl, 1./3600);



    """
    import textwrap

    print('queryname:', queryname)
    query = None

    # basic count test
    query0 = """
        SELECT
               ra
        FROM
               sdssdr7.phototag
        LIMIT 10;
        """

    query0 = """
    SELECT
       ra
    FROM
        sdssdr7.phototag
    LIMIT 10;
    """

    query0a = """
    SELECT
        COUNT(*)
    FROM
        allwise.main AS w
    LIMIT 10;
    """

    query0b = """
    SELECT COUNT(*) FROM allwise.main AS w LIMIT 10;
    """

    query0c = """
    SELECT
        COUNT(*)
    FROM
        allwise.main AS w
    WHERE
        (w.ra > 1.0 OR w.ra < 2.0) AND
        (w.dec > -1.0 AND w.dec < 1.0)
    LIMIT 10;
    """

    # join example
    query1a = """WITH x AS
    (SELECT
        ra,dec, mag_u, mag_g,mag_r,mag_i,mag_z
    from
        vst_1407.atlas_main
    where
        classification_r = -1 or
        classification_g = -1 )
    select
         mag_u,mag_g,mag_r,mag_i,mag_z,
         psfmag_u,psfmag_g,psfmag_r,psfmag_i,psfmag_z
    from
        x, lateral
    (select
        psfmag_u,psfmag_g,psfmag_r,psfmag_i,psfmag_z
    from
        sdssdr9.phototag as s
    where
        q3c_join(x.ra,x.dec,s.ra,s.dec,1./3600)
    and
        mode=1
    order by q3c_dist(x.ra, x.dec, s.ra, s.dec) limit 1 )  as y ;
    """

    # first attempt at WISE xmatch with PS1
    # This is too slow; See sub query cases.
    query1b = """
    SELECT
        w.ra, w.dec,
        w.w1mpro, w.w1sigmpro, w.w1snr,
        w.w1mag,  w.w1sigm, w.w1flg, w.w1mcor,
        w.w1mag_4, w.w1sigm_4, w.w1flg_4,
        w.w2mpro, w.w2sigmpro, w.w2snr,
        w.w2mag, w.w2sigm, w.w2flg, w.w2mcor,
        w.w2mag_4, w.w2sigm_4, w.w2flg_4,
        p.objid, p.ra, p.dec,
        p.gra, p.gdec, p.gpsfmag, p.gpsfmagerr, p.gkronmag, p.gkronmagerr,
        p.rra, p.rdec, p.rpsfmag, p.rpsfmagerr, p.rkronmag, p.rkronmagerr,
        p.ira, p.idec, p.ipsfmag, p.ipsfmagerr, p.ikronmag, p.ikronmagerr,
        p.zra, p.zdec, p.zpsfmag, p.zpsfmagerr, p.zkronmag, p.zkronmagerr,
        p.yra, p.ydec, p.ypsfmag, p.ypsfmagerr, p.ykronmag, p.ykronmagerr
    FROM
        allwise.main AS w,
        panstarrs_dr1.stackobjectthin AS p
    WHERE
        (w.w1mpro - w.w2mpro ) > 0.7 AND
        (p.zpsfmag - w.w1mpro) > 5.0 AND
        (p.ra > 1.0 OR p.ra < 2.0) AND
        (p.dec > -1.0 AND p.dec < 1.0) AND
        (w.ra > 1.0 OR w.ra < 2.0) AND
        (w.dec > -1.0 AND w.dec < 1.0) AND
        q3c_join(p.ra, p.dec, w.ra, w.dec, 3.0/3600);
    """

    # another example from Sergey
    query1c = """
    WITH a AS (
                SELECT * from table1 where column_from_table1>10
                ),
     axb AS (
                SELECT * from a, table2 AS b WHERE
                      q3c_join(a.ra, a.dec, b.ra, b.dec, 0.0002)
                )
      SELECT * FROM axb WHERE column_from_table2 < 44

    """

    # WISE region search
    query1d = """
    SELECT
        COUNT(*)
    FROM
        allwise.main AS w
    WHERE
        (w.ra > 1.0 AND w.ra < 2.0) AND
        (w.dec > -1.0 AND w.dec < 1.0);
    """


    # VST ATLAS spatial join with SDSSDR9
    query2a = """
    with x as
    (SELECT
        RA,DEC, MAG_u, mag_g,mag_r,mag_i,mag_z
    from
        vst_1407.atlas_main
    where
        classification_r = -1 or
        classification_g = -1 )
    SELECT
         mag_u,mag_g,mag_r,mag_i,mag_z,
         psfmag_u,psfmag_g,psfmag_r,psfmag_i,psfmag_z
    from
        x, lateral
    (select
        psfmag_u,psfmag_g,psfmag_r,psfmag_i,psfmag_z
    from
        sdssdr9.phototag as s
    where
        q3c_join(x.ra,x.dec,s.ra,s.dec,1./3600)
    and
        mode=1
    order by q3c_dist(x.ra, x.dec, s.ra, s.dec) limit 1 )  as y ;
    """


    # 4MOST IWG1 test region

    ra_centre = 35.5
    dec_centre = -5.5
    ra_width = 3.0
    dec_width = 3.0

    # 4MOST IW1 PS1DR1
    query3a = """
    SELECT
       *
    FROM
       panstarrs_dr1.stackobjectthin
    WHERE
        (ra BETWEEN (35.5 - 1.50) AND  (35.5 + 1.50)) AND
        (dec BETWEEN (-5.5 - 1.50) AND (-5.5 + 1.50));
    """

    query3b = """
    SELECT
       *
    FROM
       panstarrs_dr1.stackobjectthin
    WHERE
        (ra BETWEEN (%s - %s ) AND  (%s + %s )) AND
        (dec BETWEEN (%s - %s ) AND (%s + %s ));
    """ % (ra_centre, ra_width/2.0, ra_centre, ra_width/2.0,
           dec_centre, dec_width/2.0, dec_centre, dec_width/2.0)



    # WISE xmatch PS1DR1
    query3b = """
    WITH x AS (
               SELECT
                   *
               FROM
                   allwise.main
               WHERE
                   (ra > 35.5 - 1.5 AND ra < 35.5 + 1.5) AND
                   (dec > -5.5 - 1.5 AND dec < -5.5 + 1.5)
              ),
     y AS (
               SELECT
                   *
               FROM
                   x, panstarrs_dr1.stackobjectthin as z
               WHERE
                   q3c_join(x.ra, x.dec, z.ra, z.dec, 3/3600.0)
                )
      SELECT * FROM y
    """

    query3c = """
    /* WISE AGN W1W2 criterion with PS1 join and zdrop criterion */
    WITH a AS (
               SELECT
                   * from allwise.main
               WHERE
                   (w1mpro - w2mpro ) > 0.7 AND
                   (ra > 1.0 OR ra < 2.0) AND
                   (dec > -1.0 AND dec < 1.0)
              ),
     axb AS (
               SELECT
                   a.*, b.zpsfmag
               FROM
                   a,  panstarrs_dr1.stackobjectthin AS b
               WHERE
                   q3c_join(a.ra, a.dec, b.ra, b.dec, 3/3600.0)
                )
    SELECT * FROM axb WHERE (axb.zpsfmag - axb.w1mpro) > 5.0
    """


    query2b = """
    WITH a AS (
               SELECT
                   *
               FROM
                   allwise.main
               WHERE
                   (w1mpro - w2mpro ) > 0.7 AND
                   (ra > 1.0 - 1.0 AND ra < 1.0 + 1.0) AND
                   (dec > 0.0 - 1.0 AND dec < 0.0 + 1.0)
              ),
    b AS (
               SELECT
                   *
               FROM
                   panstarrs_dr1.stackobjectthin AS b
               WHERE
                   (ra > 1.0 - 1.0 AND ra < 1.0 + 1.0) AND
                   (dec > 0.0 - 1.0 AND dec < 0.0 + 1.0)
                ),
    axb AS (
               SELECT
                   a.*, b.zpsfmag
               FROM
                   a,  b
               WHERE
                   q3c_join(a.ra, a.dec, b.ra, b.dec, 3/3600.0)
                )
    SELECT * FROM axb WHERE (axb.zpsfmag - axb.w1mpro) > 5.0
    """


    # WISE with AGN selection criterion xmatch PS1DR1 and HZQ criterion
    # 3 x 3 degree region
    query2d = """
    WITH x AS (
               SELECT
                   *
               FROM
                   allwise.main
               WHERE
                   (ra > 35.5 - 1.5 AND ra < 35.5 + 1.5) AND
                   (dec > -5.5 - 1.5 AND dec < -5.5 + 1.5)
              ),
     y AS (
               SELECT
                   *
               FROM
                   x, panstarrs_dr1.stackobjectthin as z
               WHERE
                   q3c_join(x.ra, x.dec, z.ra, z.dec, 3/3600.0)
                )
    SELECT
        *
    FROM
        y
    WHERE
        (y.zpsfmag - y.w1mpro) > 2.5 AND
        (y.w1mpro - y.w2mpro ) > 0.7

    """


    # 5 x 5 degree region
    query2e = """
    WITH x AS (
               SELECT
                   *
               FROM
                   allwise.main
               WHERE
                   (ra > 35.5 - 2.5 AND ra < 35.5 + 2.5) AND
                   (dec > -5.5 - 2.5 AND dec < -5.5 + 2.5)
              ),
     y AS (
               SELECT
                   *
               FROM
                   x, panstarrs_dr1.stackobjectthin as z
               WHERE
                   q3c_join(x.ra, x.dec, z.ra, z.dec, 3/3600.0)
                )
    SELECT
        *
    FROM
        y
    WHERE
        (y.zpsfmag - y.w1mpro) > 3.5 AND
        (y.w1mpro - y.w2mpro ) > 0.7

    """



    # This is slow due to the second RA, Dec sub selection which is not needed
    query5c = """
    WITH a AS (
               SELECT
                   *
               FROM
                   allwise.main
               WHERE
                   (ra > 35.5 - 0.1 AND ra < 35.5 + 0.1) AND
                   (dec > -5.5 - 0.1 AND dec < -5.5 + 0.1)
              ),
     b AS (
               SELECT
                   *
               FROM
                   panstarrs_dr1.stackobjectthin AS b
               WHERE
                   (ra > 35.5 - 0.1 AND ra < 35.5 + 0.1) AND
                   (dec > -5.5 - 0.1 AND dec < -5.5 + 0.1)
                ),
     axb AS (
               SELECT
                   a.*, b.*
               FROM
                   a,  b
               WHERE
                   q3c_join(a.ra, a.dec, b.ra, b.dec, 3/3600.0)
                )
    SELECT * FROM axb
    """

    # query = queryname

    query = eval(queryname)

    print('query:', query)

    query = textwrap.dedent(query)

    print(query)

    return query






def getargs():
    """

    """
    import argparse

    description = '''Query WSDB database based on CASU querydqc '''
    epilog = " "
    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    default_outfile = '/tmp/tmp.fits'
    default_sqlfile = None
    default_sqlquery = None
    default_storedquery ='query0'

    parser.add_argument("--outfile",
                  default=default_outfile,
                  help="Output filename")

    parser.add_argument("--sqlfile",
                  default=default_sqlfile,
                  help="Input SQL query filename")

    parser.add_argument("--storedquery",
                  default=default_storedquery,
                  help="Stored SQL query")

    parser.add_argument("--sqlquery",
                  default=default_sqlquery,
                  help="Not implement yet;Input SQL query on command line e.g. SELECT COUNT(*) FROM table")

    parser.add_argument("--verbose",
                  action='store_true', default=False,
                  help="verbose mode")

    parser.add_argument("--debug",
                  action='store_true', default=False,
                  help="debug mode (more verbose than verbose mode)")

    args = parser.parse_args()

    return args


def clean_comments():
    """
    see https://github.com/andialbrecht/sqlparse/issues/410

    """
    import sqlparse
    # help(sqlparse.sql)
    new_query = []

    print(len(query))
    print(sqlparse.format(query, strip_comments=True))
    query = sqlparse.format(query, strip_comments=True)
    print(len(query))
    print(sqlparse.format(query, strip_comments=True))
    print(len(query))
    query = sqlparse.format(query, strip_comments=True)
    print(sqlparse.format(query, strip_comments=True))
    print(len(query))
    sys.exit()

    for statement in sqlparse.parse(query):
        # print(statement)
        new_tokens = [stm for stm in statement.tokens
                      if not isinstance(stm, sqlparse.sql.Comment)]

        print(len(new_tokens))
        for token in new_tokens:
            print(token)

        new_statement = sqlparse.sql.TokenList(new_tokens)
        new_query.append(new_statement)

    print(new_query)
    print(len(new_query))
    print(sqlparse.format("\n".join(new_query), strip_comments=True))
    sys.exit()


if __name__ == "__main__":

    print('Elapsed time(secs): ', time.time() - t0)

    # atpy support posgresl table io: http://atpy.readthedocs.io/en/stable/
    USING_atpy = False


    args = getargs()
    # help(args)
    # print('args:', args)
    print('list the args and values')
    for arg in vars(args):
        print(arg, '=', getattr(args, arg))
    print()

    # Open and read the file as a single buffer
    debug = args.debug
    sqlquery = args.sqlquery
    sqlfile = args.sqlfile
    storedquery = args.storedquery
    outfile = args.outfile

    print_table = False

    if sqlfile is not None:
        print('Reading: ',sqlfile)
        fh = open(sqlfile, 'r')
        sqlquery = fh.read()
        fh.close()
        if debug:
            print('sqlquery:', sqlquery)

    # sys.exit()

    db, host, user, password, db_table = rd_config_wsdb(table='wise')
    print('host:', host)
    print('db:', db)
    print('db_table:', db_table)
    print('user:', user)
    print('password:', password)

    if sqlquery is None:
        query = get_storedquery(queryname=storedquery)

    if sqlquery is not None:
        query = sqlquery

    if debug:
        print(query)
        for iline, line in enumerate(query.splitlines()):
            print(iline, ':', line)



    print('host:', host)
    print('query', len(query))
    print(sqlparse.format(query, strip_comments=True))
    print('getting data...')

    t1=time.time()
    result = sqlutil.get(query, db=db, host=host,
                         user=user, password=password, asDict=True)
    querytime = time.time() - t1

    if USING_atpy:
        result = atpy.Table('postgres', query=query, user=user, database='wsdb',
                            password=password, host='cappc127')

    if debug:
        help(result)

    # convert results into astropy table
    table = Table(result)
    print('len(table):', len(table))
    print('len(table.colnames):', len(table.colnames))
    print('row[0]', table[0])

    if print_table:
        table.more()

    if debug:
        help(result)

    if USING_atpy:
        result.describe()

    table.info('stats')

    print('Elapsed time(secs): ', time.time() - t0)
    print('Data downloaded:', len(table), 'rows')

    outpath = '/data/desardata4/PS/'
    outpath = '/data/4most/IWG1/'

    #outfilename = 'WISE_xmatch_PS1DR1_HZQv1.fits'
    #outfilename = 'WISE_xmatch_PS1DR1_AGN.fits'

    outfilename = 'WISE_xmatch_PS1DR1.fits'

    #outfilename = 'PS1DR1.fits'
    #outfile = outpath + outfilename


    table = table_add_sqlquery(table=table, query=query)

    table.meta['QTIME'] = querytime

    # print(table.meta)
    print('Metadata')
    for key, value in table.meta.items():
        print('{0} = {1}'.format(key, value))

    print('Save as:', outfile)
    print('Number of row:', len(table))
    print('Number of columns:', len(table.colnames))

    table.info()
    table.info('stats')
    table_metadata(table=table, verbose=verbose)

    print('Saving:', outfile)
    table.write(outfile, overwrite=True)
    print('Elapsed time(secs): ', time.time() - t0)
