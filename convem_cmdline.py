#!/usr/bin/python
import os,time,Queue,threading,thread,random
from string import *
from math import *
# End headers

# ==========================================
# Enter number of starting conditions to try
numbertries = 99
# ==========================================

endthreads = 0

class convemthread(threading.Thread):
    def __init__(self,threadID,name,q):
        self.threadID = threadID
        self.q = q
        threading.Thread.__init__(self)
    def run(self):
        CNVEM(self.q)
        
def CNVEM(q):
    while not endthreads:
        queueLock.acquire()
        if not inputqueue.empty():
            inputdata = q.get()
            inputfreqs = inputdata[0]
            cnvcount = inputdata[1]
            queueLock.release()
            
        
            for i in range(len(cnvcount)):
                cnvcount[i] = int(cnvcount[i])
            endtheloop = 0
            while endtheloop == 0:
                        numberbins = len(cnvcount)
                        maximumcn = numberbins # Expect user to enter all bins from 0 to j
                        totalindividuals = 0
                        for bins in cnvcount:
                                totalindividuals += bins
                        precision = 0.000001
                        iterations = 10000
                        pold = []
                        pnew = []
                        x = 100
                        prex = 0
                        counter = 1
                        estbinjs = []
                        haplomatrix = []
                        for eachi in range(maximumcn):
                                thislist = []
                                for eachj in range(maximumcn):
                                        thislist.append(0)
                                haplomatrix.append(thislist)
                        startingvalues = inputfreqs
                        enterowntotal = 0.0
                        for eachvalue in startingvalues:
                            enterowntotal += float(eachvalue)
                        for w in range(0,maximumcn):
                            try:
                                pold.append(float(startingvalues[w])/enterowntotal)
                                pnew.append(float(startingvalues[w])/enterowntotal)
                            except:
                                pold.append(0.0)
                                pnew.append(0.0)                        

                        # ======================================================
                        # BEGIN ITERATION
                        while abs(x-prex) > precision and counter <= iterations:
                                estbinjs = []
                                for i in range(0,maximumcn):
                                        thisbintotal = 0
                                        if i%2 == 0:
                                                # This is for even bins
                                                for j in range(0,i/2):
                                                        thisbintotal += 2.0*pnew[j]*pnew[i-j]
                                                thisbintotal += pnew[i/2]**2
                                                thisbintotal = thisbintotal * float(totalindividuals) * 2.0
                                        else:
                                                # This is for odd bins
                                                for j in range(0,int((i/2.0)-0.5)+1):
                                                        thisbintotal += 2.0*pnew[j]*pnew[i-j]
                                                thisbintotal = thisbintotal * float(totalindividuals) * 2.0
                                        estbinjs.append(thisbintotal)
                                for i in range(0,maximumcn):
                                        thisbintotal = estbinjs[i]
                                        if i%2 == 0:
                                                # This is for even bins
                                                for j in range(0,i/2):
                                                    try:
                                                        haplomatrix[j][i-j] = (float(cnvcount[i])/(totalindividuals))*((2.0*pnew[j]*pnew[i-j])/(thisbintotal/totalindividuals*2.0))
                                                    except ZeroDivisionError:
                                                        haplomatrix[j][i-j] = 0.0
                                                try:
                                                    haplomatrix[i/2][i/2] = (float(cnvcount[i])/(totalindividuals))*((pnew[i/2]**2)/(thisbintotal/totalindividuals*2.0))
                                                except ZeroDivisionError:
                                                    haplomatrix[i/2][i/2] = 0.0
                                        else:
                                                # This is for odd bins
                                                for j in range(0,int(i/2.0-0.5)+1):
                                                    try:
                                                        haplomatrix[j][i-j] = (float(cnvcount[i])/(totalindividuals))*((2.0*pnew[j]*pnew[i-j])/(thisbintotal/totalindividuals*2.0))
                                                    except ZeroDivisionError:
                                                        haplomatrix[i][i-j] = 0.0

                                for eachline in haplomatrix:
                                        thisline = []
                                        for eachcell in eachline:
                                                thisline.append("%0.4f" % eachcell)

                                # Now we've set all the expected bin contents according to est allele freqs
                                for i in range(len(pnew)):
                                        pold[i] = pnew[i]

                                for allelebin in range(0,maximumcn):
                                        sum = 0.0
                                        for j in range(0,maximumcn):
                                                for i in range(0,maximumcn):
                                                        if i == allelebin and j == allelebin:
                                                                sum += 4.0 * haplomatrix[i][j]
                                                        elif i == allelebin and j != allelebin:
                                                                sum += 2.0 * haplomatrix[i][j]
                                                        if i != allelebin and j == allelebin:
                                                                sum += 2.0 * haplomatrix[i][j]
                                        pnew[allelebin] = sum

                                prex = x
                                x = 0.0

                                for i in range(len(pold)):
                                        x += abs(pnew[i]-pold[i])

                                counter += 1
                        thisnumiterations = float(counter)
                        result = ""

                        result += """<table><tr><td><h3 style="background-color: #ddddff;  border-top: 1px solid blue;">Expected bin counts</h3><table style="text-align: center; width: 100%;" border="1" cellpadding="0" cellspacing="0"><tr><td><b>Bin</b></td><td colspan=2><b>Expected (rounded)</b></td><td colspan=2><b>Observed</b></td></tr>"""
                        
                        thischisq = 0.0
                        for i in range(len(estbinjs)):
                            expbincont = round(estbinjs[i],0)
                            obsbincont = round(cnvcount[i]*2,0)
                            if expbincont > 5 or obsbincont > 5:
                                thischisq += (float(obsbincont)-float(expbincont))**2.0/float(expbincont)
                        estimatedbincounts = []
                        for i in range(len(estbinjs)):
                                if int(round((float(estbinjs[i])/2.0))) == 0:
                                        result += "<tr><td>"+ str(i) + "</td><td width=120>" + str(int(round((float(estbinjs[i])/2.0)))) + "</td>" + "<td align='left' width=120>&nbsp;</td>"
                                        estimatedbincounts.append(str(int(round((float(estbinjs[i])/2.0)))))
                                else:
                                        result += "<tr><td>"+ str(i) + "</td><td width=120>" + str(int(round((float(estbinjs[i])/2.0)))) + "</td>" + "<td align='left' width=120><img src='../redbar.jpg' height=10 width=" + str(int((float(estbinjs[i])/2.0)/float(totalindividuals)*120)) + "></td>"
                                        estimatedbincounts.append(str(int(round((float(estbinjs[i])/2.0)))))

                                if int(str(cnvcount[i])) == 0:
                                        result += "<td width=120>" + str(cnvcount[i]) + "</td>" + "<td align='left' width=120>&nbsp;</td>" + "</tr>"
                                else:
                                        result += "<td width=120>" + str(cnvcount[i]) + "</td>" + "<td align='left' width=120><img src='../bluebar.jpg' height=10 width=" + str(int((float(cnvcount[i])/float(totalindividuals)*120))) + "></td>" + "</tr>"
                        result += "</table>"
                        result += """</td><td><h3 style="background-color: #ddddff; border-top: 1px solid blue;">Estimated allele frequencies</h3><table style="text-align: center; width: 100%;" border="1" cellpadding="0" cellspacing="0"><tr><td style='color: green; font-style:italic;'><b>Allele</b></td><td colspan=2><b>Estimated frequency (4dp)</b></td></tr>"""
                        estimatedallelefreqs = []
                        for i in range(len(pnew)):
                                if (round(pnew[i],4)) < 0.0001: #if any are missing or showing that shouldn't be, this 0.0001 needs changing up or down 1dp
                                        result += "<tr><td style='color: green; font-style:italic;'>"+ str(i) + "</td><td>&nbsp;" + str(round(pnew[i],4)) + "&nbsp;</td><td align='left' width=300>&nbsp;</td></tr>"
                                        estimatedallelefreqs.append(str(round(pnew[i],4)))
                                else:
                                        result += "<tr><td style='color: green; font-style:italic;'>"+ str(i) + "</td><td>&nbsp;" + str(round(pnew[i],4)) + "&nbsp;</td><td align='left' width=300><img src='../greenbar.jpg' height=10 width=" + str(int(round(pnew[i],4)*300)) + "></td></tr>"
                                        estimatedallelefreqs.append(str(round(pnew[i],4)))
                        result += "</table>"

                        result += """</td></tr></table><h3 style="background-color: #ddddff; border-top: 1px solid blue;">Estimated genotype frequencies</h3><table style="text-align: center; font-size: 7pt; width: 100%;" border="1" cellpadding="0" cellspacing="0">"""
                        firstrow = haplomatrix[0]
                        result += """<tr><td>&nbsp;</td>"""
                        for i in range(len(firstrow)):
                                result += """<td>""" + str(i) + """</td>"""
                        result += """</tr>"""
                        n = 0
                        for eachline in haplomatrix:
                                result += """<tr><td>""" + str(n) + """</td>"""
                                i = 0
                                for eachcell in eachline:
                                        if i >= n:
                                                realvalue = eachcell*4.0
                                                if round(realvalue,1) >= 0.1:
                                                    result += """<td bgcolor=#ffff77>%0.3f</td>""" % realvalue
                                                elif round(realvalue,2) >= 0.01:
                                                    result += """<td bgcolor=#ffffaa>%0.3f</td>""" % realvalue
                                                elif round(realvalue,3) >= 0.001:
                                                    result += """<td bgcolor=#ffffcc>%0.3f</td>""" % realvalue
                                                else:
                                                    result += """<td>%0.3f</td>""" % realvalue
                                        else:
                                                result += """<td>&nbsp;</td>"""
                                        i += 1
                                result += """</tr>"""
                                n += 1
                        result += "</table>"

                        # New unformatted output
                        observedcopynumbers = []
                        bins = []
                        n = 0
                        for item in cnvcount:
                            observedcopynumbers.append(str(item))
                            bins.append(str(n))
                            n+=1
                        result = ""
                        result+= "Copy number" + "\t" + "\t".join(bins) + "\n"
                        result+= "Observed count" + "\t" + "\t".join(observedcopynumbers) + "\n"
                        result+= "Expected count" + "\t" + "\t".join(estimatedbincounts) + "\n"
                        result+= "Allele freq" + "\t" + "\t".join(estimatedallelefreqs) + "\n\n"

                        sumest = 0.0
                        sumcount = 0.0
                        for value in estbinjs:
                            sumest += value
                        sumest = sumest/2.0
                        for value in cnvcount:
                            sumcount += float(value)
                        if sumest < sumcount-1:
                            cnvcount.append(0)
                        else:
                            endtheloop = 1
            outputqueue.put([result,pnew,thischisq,thisnumiterations])

        else:
            queueLock.release()

def pdifference(pnew,ptest):
    best = 100.0
    bestnum = 0
    try:
        if len(ptest) == 0:
            #print "short"
            return [1,1.0]
        else:
            for testcount in range(len(ptest)):
                x = 0.0
                for i in range(max(len(ptest[testcount]),len(pnew))):
                    try:
                        testx = abs(round(pnew[i],2)-round(ptest[testcount][i],2))
                        if testx > x:
                            x = testx
                    except:
                        pass
                if x < best:
                    best = x
                    bestnum = testcount
            return [bestnum,best]
    except:
        return [1,1.0]
            
try:
    inputdata = raw_input("Enter CNV bin counts as comma-delimited list, starting with bin zero:\n")
    cnvcount = inputdata.split(",")
    inputqueue = Queue.Queue(numbertries+1)
    outputqueue = Queue.Queue(numbertries+1)
    queueLock = threading.Lock()
    threads = []
    freqset = []
    queueLock.acquire()
    for countb in range(len(cnvcount)):
        freqset.append(1.0/float(len(cnvcount)))
    inputqueue.put([freqset,cnvcount])
    for counta in range(numbertries):
        freqset = []
        rawrandomnos = []
        total = 0.0

        for countb in range(len(cnvcount)):
            x = 0.0
            while x==0.0:
                x = random.random()
            rawrandomnos.append(x)
            total += x
        for item in rawrandomnos:
            freqset.append(item/total)
        inputqueue.put([freqset,cnvcount])

    queueLock.release()

    threadID = 1
    for countc in range(numbertries+1):
        thread = convemthread(threadID,str(threadID),inputqueue)
        thread.start()
        threads.append(thread)
        threadID += 1

    while not inputqueue.empty():
        pass
    endthreads = 1
    for thread in threads:
        thread.join()
        
    resultset = []
    result = []
    while not outputqueue.empty():
        resultset.append(outputqueue.get())
    ptest = []    
    countn = 1
    resultcount = {"1":0}
    resultchi = {}
    resultiterations = {}
    for item in resultset:
        testdiff = pdifference(item[1],ptest)
        if testdiff[1]>0.1:
            result.append("""\nResult """ + str(countn) + """\n"""+ str(item[0]) + "\n===========================================\n")
            ptest.append(item[1])
            resultchi[str(countn)] =  str(round(item[2],2))
            resultiterations[str(countn)] = [item[3]]
            resultcount[str(countn)] = 1
            countn += 1
        else:
            resultcount[str(testdiff[0]+1)] += 1
            resultiterations[str(testdiff[0]+1)].append(item[3])
            
    print """Possible solutions:\n"""
    iterationset = resultiterations.keys()
    for item in iterationset:
        itlist = resultiterations[item]
        counter = 0.0
        total = 0.0
        for value in itlist:
            total += value
            counter += 1.0
        resultiterations[item] = str(int(round(total/counter,0)))
    
    
    possibleset = resultcount.keys()
    possibleset.sort()
    for item in possibleset:
        print "\nResult",str(item),": ",resultcount[item],"starts (chi^2 =",resultchi[item]+", mean",resultiterations[item],"iterations)","\n"
    print """\n=======================\n"""
    for item in result:
        print item


except:
    print "Error"
