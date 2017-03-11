"""

Requires:
https://github.com/segasai/astrolibpy/tree/master/utils


this is the easyaccess git commit that adds the SQL to the FITS format
output file:
https://github.com/mgckind/easyaccess/commit/aa1ea0adfc059993ff5dfb72fd156237d5088b4d

fits[1].write_history('Created by easyaccess ' + version.__version__ + ' on ' + created)
fits[1].write_comment('Query = ' + query)

eg

HISTORY Created by easyaccess 1.3.1 on 2016-Dec-17 02:18:55                     COMMENT Query = SELECT *  FROM     SVA1_COADD WHERE     (ra BETWEEN 33.0 AND  38COMMENT .0) AND     (dec BETWEEN -8.0 AND -3.0)

better if something like

HISTORY Created by easyaccess 1.3.1 on 2016-Dec-17 02:18:55                     COMMENT STARTQUERY
COMMENT SELECT
COMMENT     *
COMMENT FROM
COMMENT     SVA1_COADD
COMMENT WHERE
COMMENT     (ra BETWEEN 33.0 AND  38.0) AND
COMMENT     (dec BETWEEN -8.0 AND -3.0)
COMMENT ENDQUERY


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
import os
import sys
import time
import datetime
import argparse

# 3rd party functions
import atpy
from astropy.table import Table
from astropy.io import fits
import pandas as pd

# private functions
sys.path.append('/home/rgm/soft/python/lib/')
from librgm import sqlutil
# help(sqlutil)
from librgm.wsdb import rd_config_wsdb
from librgm.table_metadata import table_metadata

print(atpy.__version__)
# print(sqlutil.__version__)

t0=time.time()
# print time to load all the python functions
print('Elapsed time(secs): ', time.time() - t0)


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
    SELECT
        COUNT(*)
    FROM
        allwise.main AS w
    LIMIT 10;
    """

    query0c = """
    SELECT
        COUNT(*)
    FROM
        allwise.main AS w
    WHERE
        (w.ra > 1.0 AND w.ra < 2.0) AND
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


def cutquery(query, length):
    """
    Return query in a list of fixed sized character strings

    suggest length = 70 so can be stored in FITS format COMMENT
    from easyaccess codebase
    """
    return [query[0 + i:length + i] for i in range(0, len(query), length)]


def table_metadata_sql(table, sqlfile=None, sqlquery=None):
    """
    Add SQL query information into Astropy table metadata and
    can be output in a FITS header


    """
    # various ways to store the SQL query information in the FITS header
    if sqlfile is not None:
        table.meta['SQLFILE'] = sqlfile

    # table.meta['SQL'] = query
    # create list to store sql query
    table.meta['SQLSTART'] = 'START:SQLQUERY'
    isql = -1
    for iline, line in enumerate(query.splitlines()):
        isql = isql + 1
        print(iline, ':', line)
        # cat and wrap long sql lines at 67 chars since
        # FITS has a max total record length of 80 and 13 characters
        # are needed for the required information
        # needs 12345678 = ''
        #       1234567890123
        cutresult = cutquery(line, 67)
        print('cutquery:', len(cutresult))
        for icut, record in enumerate(cutresult):
            isql = isql + icut
            if isql <= 99998:
                table.meta['SQL' + '{:d}'.format(isql+1)] = record
            if isql > 99998:
                table.meta['SQL' + '{:d}'.format(99999)] = record

    table.meta['SQLEND'] = 'END:SQLQUERY'
    table.meta['QTIME'] = querytime

    return table


def get_args():
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


if __name__ == "__main__":

    outpath = '/data/desardata4/PS/'
    outpath = '/data/4most/IWG1/'

    #outfilename = 'WISE_xmatch_PS1DR1_HZQv1.fits'
    #outfilename = 'WISE_xmatch_PS1DR1_AGN.fits'

    outfilename = 'WISE_xmatch_PS1DR1.fits'

    #outfilename = 'PS1DR1.fits'
    #outfile = outpath + outfilename



    print('Elapsed time(secs): ', time.time() - t0)

    # atpy support posgresl table io: http://atpy.readthedocs.io/en/stable/
    USING_atpy = False
    debug = False

    args = get_args()
    # help(args)
    # print('args:', args)
    print('list the args and values')
    for arg in vars(args):
        print(arg, '=', getattr(args, arg))
    print()

    # Open and read the file as a single buffer
    sqlquery = args.sqlquery
    sqlfile = args.sqlfile
    storedquery = args.storedquery
    outfile = args.outfile

    print('Output file:', outfile)

    outpath = os.path.dirname(outfile)
    print('Output file path:', outpath)

    if not os.path.exists(outpath):
       print('Output directory does not exist:', outpath)
       print('Creating:', outpath)
       os.makedirs(outpath)

    # key=raw_input("Enter any key to continue: ")

    if sqlfile is not None:
        print('Reading: ',sqlfile)
        fh = open(sqlfile, 'r')
        sqlquery = fh.read()
        fh.close()
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

    print(query)

    for iline, line in enumerate(query.splitlines()):
        print(iline, ':', line)

    print('host:', host)
    print('query set')
    print('getting data...')

    t1 = time.time()
    result = sqlutil.get(query, db=db, host=host, strLength=64,
                         user=user, password=password, asDict=True)
    querytime = time.time() - t1
    print('Elapsed time(secs): ', time.time() - t0)

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

    print_table = False
    if print_table:
        table.more()


    if debug:
        help(result)

    if USING_atpy:
        result.describe()

    table.info('stats')

    print('Elapsed time(secs): ', time.time() - t0)
    print('Data downloaded:', len(table), 'rows')


    table_metadata_sql(table, sqlfile=sqlfile, sqlquery=sqlquery)

    print('Table metadata')
    for key, value in table.meta.items():
        print('{0} = {1}'.format(key, value))

    print('Save as:', outfile)
    print('Number of row:', len(table))
    print('Number of columns:', len(table.colnames))

    table_metadata(table=table, verbose=True)

    table.write(outfile, overwrite=True)
    print('Elapsed time(secs): ', time.time() - t0)
