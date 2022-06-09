import pandas
import os
import pyreadr

# annotationDBfiles = glob.glob('data/*.tsv')[:-1]; annotationDBfile = annotationDBfiles[0]
# organisms = [10090, 9606]
# for annotationDBfile in annotationDBfiles:
#     annotationDBdf = pandas.read_csv(annotationDBfile,sep='\t')
#     annotationDBdf = annotationDBdf.loc[annotationDBdf.organism.isin(organisms),]
#     annotationDBdf.to_csv(annotationDBfile,sep='\t',index=False)

def getDorotheaDB(dorotheaDBfile = 'data/dorothea.tsv', org = 9606, confidence='C'):
    if os.path.exists(dorotheaDBfile):
        dorotheaDBdf = pandas.read_csv(dorotheaDBfile,sep='\t')
        dorotheaDBdf = dorotheaDBdf.loc[dorotheaDBdf.org == org]
        dorotheaDBdf = dorotheaDBdf.loc[dorotheaDBdf.confidence <= confidence]
        return(dorotheaDBdf)
    else:
        dorotheaDBurls = ['https://github.com/saezlab/dorothea/raw/master/data/dorothea_mm.rda',
                          'https://github.com/saezlab/dorothea/raw/master/data/dorothea_hs.rda']
        dorotheaDBfiles = []
        for dorotheaDBurl in dorotheaDBurls:
            dorotheaDBfile = dorotheaDBurl.split('master/')[1]
            dorotheaDBfiles.append(dorotheaDBfile)
            if os.path.exists(dorotheaDBfile) == False:
                cmd = ['wget',dorotheaDBurl,'-P','data/']
        dorotheaDBfile = dorotheaDBfiles[0]
        objName = os.path.splitext(os.path.basename(dorotheaDBfile))[0]
        dorotheaDB = pyreadr.read_r(dorotheaDBfile)
        dorotheaDBdf = dorotheaDB[objName]
        keycols = ['tf','target','confidence']
        dorotheaDBdf = dorotheaDBdf[keycols]
        org = 10090 if 'mm' in objName else 9606
        dorotheaDBdf['org'] = org
        for dorotheaDBfile in dorotheaDBfiles[1:]:
            objName = os.path.splitext(os.path.basename(dorotheaDBfile))[0]
            dorotheaDB = pyreadr.read_r(dorotheaDBfile)
            org = 10090 if 'mm' in objName else 9606
            tmp_dorotheaDBdf = dorotheaDB[objName][keycols]
            tmp_dorotheaDBdf['org'] = org
            dorotheaDBdf = pandas.concat([dorotheaDBdf, tmp_dorotheaDBdf])
        dorotheaDBdf.drop_duplicates(inplace = True)
        dorotheaDBdf.to_csv('data/dorothea.tsv',sep='\t',index=False)
        return(dorotheaDBdf)
