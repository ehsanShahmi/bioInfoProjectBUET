#from vpnocchio import VPN, init_logging
from threading import Thread
import pandas as pd
from os import listdir
from os.path import isfile, join
from pprint import pprint
import requests
import datetime
import json
import glob
import ftplib
import GEOparse
import os
import re

'''
def use_ovpn(task, task_args,
             ovpn_config_dir='/home/arif/ovpn-configs/VPNBook.com-OpenVPN-US1',
             min_time_before_reconnect = 30,
             log_on = False,
             credentials=[('vpnbook', '533d2ve', 'US')]):
    if log_on:
        init_logging()

    VPN.conf_dir = ovpn_config_dir
    VPN.min_time_before_reconnect = min_time_before_reconnect

    #credentials = [('vpnbook', '533d2ve', 'US')]

    def do_something(t, ta, *args):
        vpn = VPN(*args)
        for one in range(1):            
            # it has requests inside
            response = vpn.get('http://ip.barjomet.com')
            vpn.log.info('Current IP: %s',  response.text)
            t(**ta)
            #vpn.new_ip()
            vpn.disconnect()

    for username, password, match_config_name in credentials:
        do_something(task, task_args, username, password, match_config_name)
        #Thread(target=do_something,
               #args=(t=task,ta=task_args,
                     #username,
                     #password,
                     #match_config_name)).start()


'''
def collect_gene_profiles(gene_sym, output_dir):
    db = "geoprofiles"
    term = "("+gene_sym+"[Gene Symbol]) AND Homo sapiens[Organism]"
    retmax = "1000000"
    usehistory = "y"
    retmode = "json"
    retstart = "0"
    
    request_count = 1 # if #idlist > retmax, we'll need multiple requests

    response = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=" + db + "&term="+term+"&retmax="+retmax+"&usehistory="+usehistory+"&retmode="+retmode)
    
    if(response.ok):
        with open(os.path.join(output_dir, gene_sym+ '_geoprofiles_list'+ '.json'), 'w') as outfile:
            outfile.write(response.text)
            #json.dump(response.content.decode('utf-8'), outfile)
            
        jData = json.loads(response.content.decode('utf-8'))            
        #print(json.dumps(jData, indent=4, sort_keys=True))
        #print(json.dumps(jData["esearchresult"]["idlist"], indent=4, sort_keys=True))
        print(datetime.datetime.now().strftime('%d %b %Y %H:%M:%S') + ' | DONE | List of geoprofiles downloaded for | >' + gene_sym + '<')
    else:
        response.raise_for_status()



#collect_gene_profiles('A1CF', './geo-data')
#use_ovpn(collect_gene_profiles, {'gene_sym': 'SHH', 'output_dir': './geo-data'})


def collect_gene_profiles_for_711_genes():
    already_done_genes = set()

    already_done_genes.add('SEPT5#') # problematic
    already_done_genes.add('SEPT6#')
    already_done_genes.add('SEPT9#')
    
    with open('./geo-data/geo-profiles/log') as log_file:
        log_data = log_file.readlines()
        for log in log_data:
            #print(log)
            matches = re.search(r' >(.*)<$', log)
            if matches:
                #print(matches.group(1))
                already_done_genes.add(matches.group(1))
            else:
                print ("No match!!")
        #print(log_data)

    print('Previously done: ' + str(len(already_done_genes)))
        
    cancer_711_df = pd.read_csv('./NCG6_tsgoncogene.tsv', sep='\t', encoding='utf-8')
    gene_syms_711 = cancer_711_df.iloc[:, 1]
    count = 0;
    for gene_sym in  gene_syms_711:
        if gene_sym not in already_done_genes:            
            #print(gene_sym)
            #use_ovpn(collect_gene_profiles, {'gene_sym': gene_sym, 'output_dir': './geo-data/geo-profiles'})
            collect_gene_profiles(gene_sym, './geo-data/geo-profiles')
            count += 1
    print('DONE: ' + str(count))

#collect_gene_profiles_for_711_genes()

def get_unique_gene_profiles():
    dir_path = './geo-data/geo-profiles'
    #geneprofile_files = [f for f in listdir(dir_path) if isfile(join(dir_path, f))]

    #for geneprofile_file in geneprofile_files:
        #if geneprofile_file.
        #print(geneprofile_file)

    #print(glob.glob('/mnt/drive01/MSc/BioInfo/project/geo-data/geo-profiles/*_geoprofiles_list.json'))
    geneprofiles_files = glob.glob('./geo-data/geo-profiles/*_geoprofiles_list.json')
    print('File found: ' + str(len(geneprofiles_files)))

    id_set = set()
    id_list = []
    for geneprofiles_file in geneprofiles_files:

        with open(geneprofiles_file) as geneprofile_json:
            print('Reading ' + geneprofiles_file)
            data = json.load(geneprofile_json)
            ids = data["esearchresult"]["idlist"]
            #print(ids)
            for id in ids:
                #print(id)
                id_set.add(id.encode('utf-8').strip())
                id_list.append(id.encode('utf-8').strip())

    print('Unique gene profile found: ' + str(len(id_set)) + ' from:' + str(len(id_list)))
    # Unique gene profile found: 2615121 from:2632666

    with open('./geo-data/geo-profiles/unique-profile-ids-02.json', 'w')  as out_file:
        json.dump(list(id_set), out_file)

#get_unique_gene_profiles()



def get_geneprofile_summary(profile_id, output_dir):
    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=geoprofiles&id='+profile_id+'&retmode=json')

    if(response.ok):
        profile_ids = [x.strip() for x in profile_id.split(',')]
        
        with open(os.path.join(output_dir, profile_ids[0]+ '_geoprofile_summary'+ '.json'), 'w') as outfile:
            outfile.write(response.text.encode('utf8'))
        
        for p_id in profile_ids:
            print(datetime.datetime.now().strftime('%d %b %Y %H:%M:%S') + ' | DONE | geoprofile summary downloaded for | >' + p_id + '<')
    else:
        response.raise_for_status()

#get_geneprofile_summary('132669597')



def get_geoprofile_summaries(id_json_file, output_dir):
    already_done_profiles = set()
    
    with open('./geo-data/geo-profile-summaries/log') as log_file:
        log_data = log_file.readlines()
        for log in log_data:
            #print(log)
            matches = re.search(r' >(.*)<$', log)
            if matches:
                #print(matches.group(1))
                already_done_profiles.add(matches.group(1))
            else:
                print ("No match!!")
        #print(log_data)

    print('Previously done: ' + str(len(already_done_profiles)))
    
    with open(id_json_file, 'r')  as profile_ids_file:
        ids_all = json.load(profile_ids_file)
        ids = [x for x in ids_all if x not in already_done_profiles]
        print('Previously done: ' + str(len(already_done_profiles)) + ' or' + str(len(ids_all) - len(ids))  + ' from ' + str(len(ids_all)))

        def chunks(l, n):
            """Yield successive n-sized chunks from l."""
            for i in xrange(0, len(l), n):
                yield l[i:i + n]

        splits = chunks(ids, 250)
        #print(type(splits.next()))

        #id_param = ''
        #split = splits.next()
        #print(str(len(split)))
        #c = 0;
        #for id in split:
            #id_param += id + ','
        #    c += 1
        #    print(id)
        #print('param: ' + id_param.rstrip(','))
        #print('Count: ' + str(c))
        #get_geneprofile_summary(id_param.rstrip(','), output_dir)
        
        for split in splits:
            id_param = ''
            for id in split:
                id_param += id + ','
            #print(id_param.rstrip(','))
            get_geneprofile_summary(id_param.rstrip(','), output_dir)
        
        #for id in ids:
            #get_geneprofile_summary(id, output_dir)

#get_geoprofile_summaries('./geo-data/unique-profile-ids-utf8.json', './geo-data/geo-profile-summaries')
#use_ovpn(get_geoprofile_summaries, {'id_json_file':'./geo-data/unique-profile-ids-utf8.json', 'output_dir':'./geo-data/geo-profile-summaries'})

def get_unique_gds_ids():
    #dir_path = './geo-data/geo-profile-summaries-backup'
        
    #print(glob.glob('/mnt/drive01/MSc/BioInfo/project/geo-data/geo-profiles/*_geoprofiles_list.json'))
    geneprofile_summary_files = glob.glob('./geo-data/geo-profile-summaries/*_geoprofile_summary.json')
    print('File found: ' + str(len(geneprofile_summary_files)))

    id_set = set()
    id_list = []
    for geneprofile_summary_file in geneprofile_summary_files:

        with open(geneprofile_summary_file) as geneprofile_summary_json:
            print('Reading ' + str(geneprofile_summary_json))
            data = json.load(geneprofile_summary_json)
            result = data["result"]
            geneprofile_uids = result["uids"]
            #print(ids)
            for geneprofile_uid in geneprofile_uids:
                gds_id = result[str(geneprofile_uid)]['gds']
                #print('GDS: ' + gds_id)
                id_set.add(gds_id.encode('utf-8').strip())
                id_list.append(gds_id.encode('utf-8').strip())

    print('Unique GDS found: ' + str(len(id_set)) + ' from: ' + str(len(id_list)))
    # Unique GDS found: 1772 from: 1412305
    # Unique GDS found: 1772 from: 1711805

    with open('./geo-data/unique-gds-ids.json', 'w')  as out_file:
        json.dump(list(id_set), out_file)

#get_unique_gds_ids()


def get_gds_summary(gds_ids, output_dir):
    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id='+gds_ids+'&retmode=json')

    if(response.ok):
        gds_ids_splitted = [x.strip() for x in gds_ids.split(',')]
        
        with open(os.path.join(output_dir, gds_ids_splitted[0]+ '_gds_summary'+ '.json'), 'w') as outfile:
            outfile.write(response.text.encode('utf8'))
        
        for g_id in gds_ids_splitted:
            print(datetime.datetime.now().strftime('%d %b %Y %H:%M:%S') + ' | DONE | GDS summary downloaded for | >' + g_id + '<')
    else:
        response.raise_for_status()

#get_gds_summary('6100,6177', './geo-data/gds-summaries')


def get_gds_summaries(gds_ids_file):

    def chunks(l, n):
        """Yield successive n-sized chunks from l."""
        for i in xrange(0, len(l), n):
            yield l[i:i + n]

    with open(gds_ids_file) as gds_file:
        ids = json.load(gds_file)

        splits = chunks(ids, 250)

        for split in splits:
            #print(','.join(split))
            get_gds_summary(','.join(split), './geo-data/gds-summaries')
        
#get_gds_summaries('./geo-data/unique-gds-ids.json')        

def extract_ftp_link_for_gds(gds_summary_dir):
    gds_summary_files = glob.glob('./geo-data/gds-summaries/*_gds_summary.json')

    ftp_links = set()

    for gds_summary_file in gds_summary_files:
        with open(gds_summary_file) as gds_summary_json:
            print('Reading ' + str(gds_summary_json))
            data = json.load(gds_summary_json)
            result = data["result"]
            gds_uids = result["uids"]

            for gds_uid in gds_uids:
                ftp_link = result[str(gds_uid)]['ftplink']
                print(ftp_link)
                ftp_links.add(ftp_link.encode('utf-8').strip())

    print('#FTP links for GDS found: ' + str(len(ftp_links)))

    with open('./geo-data/gds-ftp-links.json', 'w')  as out_file:
        json.dump(list(ftp_links), out_file)

#extract_ftp_link_for_gds('./geo-data/gds-summaries')


def get_links_for_GDS_SOFT_files(ftp_links_file):

    soft_file_download_links = set()

    with open(ftp_links_file) as ftp_links_json:
        links = json.load(ftp_links_json)
        
        for link in links:            
            soft_file_download_links.add(link + 'soft/' + link.split('/')[-2]+'.soft.gz')

        print('#links generated: ' + str(len(soft_file_download_links)))

        with open('./geo-data/gds-ftp--soft-download-links.json', 'w')  as out_file:
            json.dump(list(soft_file_download_links), out_file)

#get_links_for_GDS_SOFT_files('./geo-data/gds-ftp-links.json')


def print_links(link_file):
    with open(link_file) as links_json:
        links = json.load(links_json)
        
        for link in links:
            print(link)

#print_links('./geo-data/unique-gds-ids.json')


def download_from_ftp(links_files='ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/GDS1481/soft/GDS1481.soft.gz'):
    ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login('anonymous', 'anonymous@mail.com')
    
    ftp.cwd('/geo/datasets/GDS1nnn/GDS1481/soft')

    output_dir = './temp'
    os.chdir(output_dir)

    file_name = links_files.split('/')[-1]
    print(file_name)
    file = open(file_name, 'wb')
    ftp.retrbinary('RETR ' + file_name, file.write)
    file.close()

#download_from_ftp()



def get_gds_soft_files(gds_ids_file, output_dir):
    def get_gds_with_geo_parse(gds_id, output_dir):
        geo_ass = "GDS"+gds_id
        print(geo_ass)
        gse = GEOparse.get_GEO(geo=geo_ass, destdir=output_dir)

    with open(gds_ids_file, 'r') as ids_file:
        ids = ids_file.readlines()
        for id in ids:
            get_gds_with_geo_parse(id.rstrip(), output_dir)

get_gds_soft_files('./gds_ids_aa', './gds-soft')
