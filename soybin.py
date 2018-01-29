/*
* Author: Juned Ahmed
* Description: A simple genomic data scraper from https://phytozome.jgi.doe.gov
*/


from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support.expected_conditions import staleness_of
from selenium.webdriver.support import expected_conditions as ec
from selenium import webdriver
from bs4 import BeautifulSoup as bs
from Bio.SeqUtils.ProtParam import ProteinAnalysis as pa
from collections import deque

import csv
import time
import sys
import re  
reload(sys)  
sys.setdefaultencoding('utf-8')

def parse(url):
	#wd = webdriver.Chrome()
	wd = webdriver.PhantomJS()
	wd.set_window_size(1280, 800)
	wd.get(url)
	time.sleep(3)
	seq = wd.find_element_by_xpath("//td[contains(., 'Sequences')]").click()
	time.sleep(2)
	sall = wd.find_element_by_xpath("//td[text()='Show all']").click()
	hp = wd.page_source
	#wd.save_screenshot('out.png')
	wd.quit()
	data = bs(hp, "html.parser")
	return data

def geninfo(data):
	gendata = {}
	for tr in data.find("table", "infoBlock").find_all("tr"):
		td0 = ""
		for td in tr.find_all("td"):
			if (not td0):
				td0 = td.text
				gendata[td0] = None
			else:
				if td0 == "Other transcripts":
					gendata[td0] = td.find_all("a")
				else:
					gendata[td0] = td.text 
	cdscoord = re.split("\:|\..| " ,gendata["Location:"])
	tn = gendata["Transcript Name"].split(" ")
	cdscoord.append(tn[0])
	otherts = []
	try:
		if (gendata["Other transcripts"]):
			for ot in gendata["Other transcripts"]:
				otherts.append(ot.get("href"))
			return (cdscoord, otherts)
	except:
		otherts = []
		return (cdscoord, otherts)



def seqinfo(data):
	seqdic = {}
	for td in data.find_all("td", "sectionHeaderopened"):
		dt = td.text.split('[')
		seqdic[dt[0].strip()] = int(filter(str.isdigit, str(dt[1])))
	return seqdic



def protinfo(data):
	pepinfo = {}
	pepseq = ""
	for pr in data.find("div", id='PeptideSequence').find('span').find_all('span'):
		pepseq += pr.text
	anaseq = pa(str(pepseq).strip("*"))
	pepinfo["PP"] = len(str(pepseq).strip("*"))
	pepinfo["MW"] = anaseq.molecular_weight()
	pepinfo["pI"] = anaseq.isoelectric_point()

	return pepinfo



def seqdata(data):
	pepseq = ""
	cdsseq = ""
	genseq = ""
	matchpep = re.compile('\>\w+\.\w+\.\d')
	matchcds = re.compile('\>\w+\.\w+\.\d\s(CDS)')
	matchgen = re.compile('\>\w+\.\w+\s\|\s\w+\:\d+\.+\d+\s(forward)|(reverse)$')
	for cds in data.find_all("div", class_='exonBlock'):
		if (re.search(matchcds, cds.text)):
			cdsseq += re.search(matchcds, cds.text).group() + "\n"
			for c in cds.find_all('span'):
				cdsseq += c.text
		elif (re.search(matchgen, cds.text)):
			genseq += re.search(matchgen, cds.text).group() + "\n"
			for g in cds.find_all('span'):
				genseq += g.text
	pepseq += re.search(matchpep, data.find('div', id='PeptideSequence').text).group() + ' Peptide\n'
	for pr in data.find("div", id='PeptideSequence').find('span').find_all('span'):
		pepseq += pr.text
	pepseq = pepseq.strip("*").strip("%")
	return (genseq, cdsseq, pepseq)


def datacat(dt, cnt):
	data = {"SL No": None, 
			"Gene Name": None, 
			"Locus": None, 
			"Chr. No": None,
			"Strand": None, 
			"CDS Coord (5' to 3')": None, 
			"Gene (bp)": None, 
			"CDS (bp)": None,
			"PP (aa)": None, 
			"MW (kDa)": None, 
			"pI": None}
	gi, ots = geninfo(dt)
	si = seqinfo(dt)
	pi = protinfo(dt)
	data["SL No"] = cnt + 1
	data["Gene Name"] = "SlGSTU5"
	data["Locus"] = gi[4]
	data["Chr. No"] = int(filter(str.isdigit, str(gi[0])))
	if (gi[3] == u"forward"):
		data["Strand"] = u"+"
	else:
		data["Strand"] = u"-"
	data["CDS Coord (5' to 3')"] = gi[1] + "-" + gi[2]
	data["Gene (bp)"] = si["Genomic Sequence"]
	data["CDS (bp)"] = si["CDS Sequence"]
	data["PP (aa)"] = pi["PP"]
	data["MW (kDa)"] = pi["MW"]
	data["pI"] = pi["pI"] 
	return (data, ots)

def csvwrite(data)

def main():

	locusName = []
	tmpurl = deque()
	with open(sys.argv[1], 'r') as fl:
		ln = re.findall('[A-Z]\D+\d+\w\d{6}', fl.read())
		locusName = list(set(ln))
	print locusName, len(locusName)

	burl = "https://phytozome.jgi.doe.gov/pz/portal.html#!gene?organism=Gmax&searchText=locusName:"
	counter = 0
	for l in locusName:
		try:
			url = burl+l
			dt = parse(url)
			data, ots = datacat(dt, counter)
			gd, cd, pd = seqdata(dt)
			if (ots != []):
				for ot in ots:
					counter += 1
					ourl = "https://phytozome.jgi.doe.gov/pz/portal.html"
					uri = ourl + ot
					dt = parse(uri)
					dat, ots = datacat(dt, counter)
					print "ots: ", dat
					with open("SlGSTU5.csv", "a+") as csvf:
						sn = csv.Sniffer()
						csvf.seek(0)
						csvf.seek(0)
						csvf.seek(0)
						try:
							hh = sn.has_header(csvf.read(2048))
						except:
							hh = False
						dw = csv.DictWriter(csvf, dat.keys(), dialect="excel")
						dw.writeheader()
						if (csvf.read() != "" and hh):
							pass
						else:
							dw.writeheader()
						dw.writerow(dat)
					print "ots csv write"
					gd, cd, pd = seqdata(dt)
					with open("SEQDATA-"+sys.argv[1], "a+") as sd:
						sd.write("\n>Locus: " + dat["Locus"] + "\n\n")
						sd.write(gd + "\n")
						sd.write(cd + "\n")
						sd.write(pd + "\n")
					print "ots seqdt write"
			print "dt: ", data
			with open("SlGSTU5.csv", "a+") as csvf:
				sn = csv.Sniffer()
				csvf.seek(0)
				csvf.seek(0)
				csvf.seek(0)
				try:
					hh = sn.has_header(csvf.read(2048))
				except:
					hh = False
				dw = csv.DictWriter(csvf, data.keys(), dialect="excel")
				if (csvf.read() != "" and hh):
					pass
				else:
					dw.writeheader()
				dw.writerow(data)
				print "pts csv write"
			with open("SEQDATA-"+sys.argv[1], "a+") as sd:
				sd.write("\n>Locus: " + data["Locus"] + "\n\n")
				sd.write(gd + "\n")
				sd.write(cd + "\n")
				sd.write(pd + "\n")
			print "pts seqdata write"
			counter += 1
		except Exception as e:
			print e
			tmpurl.append(l)
			continue
		finally:
			print "Success!"



	if (tmpurl):
		for l in tmpurl:
			try:
				turl = burl+l
				dt = parse(turl)
				data, tots = datacat(dt, counter)
				gd, cd, pd = seqdata(dt)
				if (tots != []):
					for ot in tots:
						counter += 1
						tourl = "https://phytozome.jgi.doe.gov/pz/portal.html"
						turi = tourl + ot
						dt = parse(turi)
						dat, tots = datacat(dt, counter)
						print "tots: ", dat
						with open("SlGSTU5.csv", "a+") as csvf:
							sn = csv.Sniffer()
							csvf.seek(0)
							csvf.seek(0)
							csvf.seek(0)
							try:
								hh = sn.has_header(csvf.read(2048))
							except:
								hh = False
							dw = csv.DictWriter(csvf, dat.keys(), dialect="excel")
							if (csvf.read() != "" and hh):
								pass
							else:
								dw.writeheader()
							dw.writerow(dat)
						print "tots csv write"
						with open("SEQDATA-"+sys.argv[1], "a+") as sd:
							sd.write("\n>Locus: " + dat["Locus"] + "\n\n")
							sd.write(gd + "\n")
							sd.write(cd + "\n")
							sd.write(pd + "\n")
						print "tots seqdata write"
				print "dt: ", data
				with open("SlGSTU5.csv", "a+") as csvf:
					sn = csv.Sniffer()
					csvf.seek(0)
					csvf.seek(0)
					csvf.seek(0)
					try:
						hh = sn.has_header(csvf.read(2048))
					except:
						hh = False
					dw = csv.DictWriter(csvf, data.keys(), dialect="excel")
					if (csvf.read() != "" and hh):
						pass
					else:
						dw.writeheader()
					dw.writerow(data)
					print "tpts csv write"
				with open("SEQDATA-"+sys.argv[1], "a+") as sd:
					sd.write("\n>Locus: " + data["Locus"] + "\n\n")
					sd.write(gd + "\n")
					sd.write(cd + "\n")
					sd.write(pd + "\n")
				print "tpts seqdata write"
				counter += 1
			except Exception as e:
				print e
				with open("err-"+sys.argv[1], 'w+') as fl:
					fl.write(l + "\n")
				continue
			finally:
				print "Success!"


	 			
if __name__ == "__main__":
	main()

