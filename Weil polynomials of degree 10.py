#!/usr/bin/env python
# coding: utf-8

# In[119]:


#Variables
R.<x> = ZZ[]
C = ComplexField(10000)
p = 2
n = 1
q = p^n
dim = 5
deg = 2 * dim

#Change p,n in variables and below in the loads to your new value of q = p^n, up to q = 3 if you want the whole file to be able to run

#change 2 in the load the below to match the new value of q, so that it is load("q-weil_polynomials_dim4.sage")
#can use any q if you want to run everything except the last parts after def CharpolysofsimpleAVofdim5frombounds(l,p,n):

load("2-weil_polynomials_dim5.sage")

#the file that has been loaded contains the characteristic polynomials of all simple abelian varieties over F_q of dimension 5
#this data has been downloaded from the LMFDB and the polynomials have been separated based on p-rank
#The LMFDB-data uses L-polynomials (which are Weil polynomials, but with the coefficients reversed)
#Hence we convert each element in the lists to Weil polynomials

def ConvertLtoWeil(LMFDBlist):
    Weilpolys = []
    for polynomials in LMFDBlist:
        d = polynomials.degree()
        newpolycoeffs = []
        for i in range(d, -1, -1):
            newpolycoeffs.append(polynomials[i])
        Weilpolys.append(R(newpolycoeffs))
    return Weilpolys

charpolydim5prank5 = ConvertLtoWeil(data_dim5_prank5)
charpolydim5prank4 = ConvertLtoWeil(data_dim5_prank4)
charpolydim5prank3 = ConvertLtoWeil(data_dim5_prank3)
charpolydim5prank2 = ConvertLtoWeil(data_dim5_prank2)
charpolydim5prank1 = ConvertLtoWeil(data_dim5_prank1)
charpolydim5prank0 = ConvertLtoWeil(data_dim5_prank0)
charpolydim5 = charpolydim5prank5 + charpolydim5prank4 + charpolydim5prank3 + charpolydim5prank2 + charpolydim5prank1 + charpolydim5prank0

#data_dim_5_prankm is the list of data of simple abelian varieties of dimension 5 over F_q with p-rank m


# In[120]:


#Built-in function in Sage to generate all q-Weil polynomials of degree deg
L = R.weil_polynomials(deg, q)


# In[121]:


#Generate list of all possible coefficients as in Theorem 9.2.3

def Weilpolysofdeg10frombounds(q):

    l = []

    bounda1lower = math.floor(- 10 * sqrt(q)) + 1
    bounda1higher = math.ceil(10 * sqrt(q))


    for a1 in range(bounda1lower, bounda1higher):
        u2part1 = 3/25 * a1^2 + 3/2 * q
    
        u3part1 = - 2/125 * a1^3 + a1 * q / 10
    
        u4part1 = - 3/625 * a1^4 + a1^2 * q / 5 + q^2
    
        u5part1 = -4/3125 * a1^5
    
        bounda2lower = math.floor(8 * sqrt(q) * abs(a1) - 35 * q) + 1
        bounda2higher = math.floor(10/3 * u2part1) + 1
    
        for a2 in range(bounda2lower, bounda2higher):
            u2 = -u2part1 + 3/10 * a2
        
            u3part2 = 3 * a1 * a2 / 50 + u3part1
        
            u4part2 = u4part1 + 3/125 * a1^2 * a2 - 3/5 * q * a2
        
            u5part2 = u5part1 + a1^3 * (15 * q + a2)/125
        
            bounda3abs = 10 * u3part2
        
            bounda3firstvar = 1/50 * (-100/3 * u2)^(3/2)
        
            bounda3lower1 = math.ceil(bounda3abs - bounda3firstvar)
            bounda3higher1 = math.floor(bounda3abs + bounda3firstvar) + 1
        
            bounda3secondabs = 6 * sqrt(q) * a2 + 50 * q * sqrt(q)
            bounda3secondvar = - 20 * q * a1
        
            bounda3lower2 = math.floor(- bounda3secondabs + bounda3secondvar) + 1
            bounda3higher2 = math.ceil(bounda3secondabs + bounda3secondvar)
        
            bounda3lower = max(bounda3lower1, bounda3lower2)
            bounda3higher = min(bounda3higher1, bounda3higher2)
        
            for a3 in range(bounda3lower, bounda3higher):
                u3 = - u3part2 + a3 / 10
            
                u4part3 = u4part2 - 2/25 * a1 * a3
            
                u5part3 = u5part2 - a1^2 * a3 / 25 + 2 * q * a3
            
                bounda4abs = - 5 * u4part3 + 10/3 * u2^2
            
                etareal = - u2^6 / 27 - 5 * u2^3 * u3^2 + 27/2 * u3^4
                etaimaginary = 27/2 * u3 * (- u3^2 - 4/27 * u2^3)^(3/2)
            
                base = C(etareal, etaimaginary)
                eta = base.nth_root(3)
                eta2 = eta * E(3)
                eta3 = eta * E(3)^2
            
                root1 = 10 * real_part(eta)
                root2 = 10 * real_part(eta2)
                root3 = 10 * real_part(eta3)
            
                bounda4lower1 = math.ceil(bounda4abs + min([root1, root2, root3]))
                bounda4higher = math.floor(bounda4abs + median([root1, root2, root3])) + 1
            
                bounda4lower2 = math.floor(4 * sqrt(q) * abs(4 * a1 * q + a3) - 9 * q * a2 - 25 * q^2) + 1
            
                bounda4lower = max(bounda4lower1, bounda4lower2)
            
                for a4 in range(bounda4lower, bounda4higher):
                    u4 = u4part3 + a4 / 5
                
                    if u3 == 0:
                        xi1 = sqrt(-u2 + sqrt(u2^2 - u4))
                        xi2 = sqrt(-u2 - sqrt(u2^2 - u4))
                        xi3 = - xi1
                        xi4 = - xi2
                    else:
                        v2 = - u2^2 / 3 - u4
                        v3 = 2/3 * u2 * u4 - 2/27 * u2^3 - 2 * u3^2
                        if v2 == 0:
                            constbase = C(- v3)
                            y = constbase.nth_root(3) - 2/3 * u2
                        
                        else:
                            delta = v3^2 + 4/27 * v2^3
                            if delta < 0:
                                constbase = C(-v3/2, sqrt(- delta)/ 2)
                            
                            else:
                                constbase = C((- v3 + sqrt(delta))/2)
                        
                            const = constbase.nth_root(3)
                            y = const - v2 / (3 * const) - 2/3 * u2
                        
                        sqrt2y = sqrt(2 * y)
                        inroot = - 4 * u2 - 2 * y
                        inrootabs = 8 * u3 / sqrt2y
                    
                        wholesqroot1 = sqrt(inroot - inrootabs)
                        wholesqroot2 = sqrt(inroot + inrootabs)
                    
                        xi1 = (sqrt(2 * y) + sqrt(-4 * u2 - 2 * y - 8 * u3 / sqrt(2 * y))) / 2
                        xi2 = (sqrt(2 * y) - sqrt(-4 * u2 - 2 * y - 8 * u3 / sqrt(2 * y))) / 2
                        xi3 = (-sqrt(2 * y) + sqrt(-4 * u2 - 2 * y + 8 * u3 / sqrt(2 * y))) / 2
                        xi4 = (-sqrt(2 * y) - sqrt(-4 * u2 - 2 * y + 8 * u3 / sqrt(2 * y))) / 2
                
                    xis = [real(xi1), real(xi2), real(xi3), real(xi4)]
                    xis.sort()
                
                    ggi4 = - xis[0]^5 - 10/3 * u2 * xis[0]^3 - 10 * u3 * xis[0]^2 - 5 * u4 * xis[0]
                    ggi3 = - xis[1]^5 - 10/3 * u2 * xis[1]^3 - 10 * u3 * xis[1]^2 - 5 * u4 * xis[1]
                    ggi2 = - xis[2]^5 - 10/3 * u2 * xis[2]^3 - 10 * u3 * xis[2]^2 - 5 * u4 * xis[2]
                    ggi1 = - xis[3]^5 - 10/3 * u2 * xis[3]^3 - 10 * u3 * xis[3]^2 - 5 * u4 * xis[3]
                
                    bounda5abs = u5part3 + a1 / 5 * (a4 - 3 * q * a2 - 5 * q^2)
                
                    bounda5lower1 = math.ceil(bounda5abs + max([ggi2, ggi4]))
                    bounda5higher1 = math.floor(bounda5abs + min([ggi1, ggi3])) + 1
                
                    bounda5secondabs = 2 * a4 * sqrt(q) + 2 * sqrt(q) * q * a2 + 2 * q^2 * sqrt(q)
                    bounda5secondvar = -2 * q * a3 - 2 * q^2 * a1
                
                    bounda5lower2 = math.floor(-bounda5secondabs + bounda5secondvar) + 1 
                    bounda5higher2 = math.ceil(bounda5secondabs + bounda5secondvar)
                
                    bounda5lower = max(bounda5lower1, bounda5lower2)
                    bounda5higher = min(bounda5higher1, bounda5higher2)
                
                    for a5 in range(bounda5lower, bounda5higher):
                        l.append(R([q^5, a1 * q^4, a2 * q^3, a3 * q^2, a4 * q, a5, a4, a3, a2, a1, 1]))
            
    return(l)

l = Weilpolysofdeg10frombounds(q)


# In[122]:


#Check for differences between lists generated by the built-in function and the one generated by Theorem 11.2.3
def comparisonlists(l,L):
    Missing = []

    for e in L:
        if e not in l:
            Missing.append(e)
            print(e.factor())

    print("number of polynomials missed by the bounds:", len(Missing))

    extra = []
    for f in l:
        if f not in L:
            extra.append(f)
            print(f.factor())

    print("number of extra polynomials found by the bounds", len(extra))

    Missingreals = []

    for pols in Missing:
        if (pols(sqrt(q)) == 0) or (pols(- sqrt(q)) == 0):
            Missingreals.append(pols)

    print("number of the missing polynomials that have a real root:", len(Missingreals))

    realsinresult = []

    for realtest in l:
        if (realtest(sqrt(q)) == 0) or (realtest(- sqrt(q)) == 0):
            realsinresult.append(realtest)

    print("number of weil polynomials with real roots satisfying the bounds:", len(realsinresult))

comparisonlists(l,L)


# In[10]:


#The polynomials missing are exactly the ones with a real root, as described in Theorem 11.2.2 (and 11.2.1 if q is a square)


# In[124]:


#check which Weil polynomials are characteristic polynomials of simple abelian varieties over F_q, according to Theorem 11.4

def CharpolysofsimpleAVofdim5frombounds(l,p,n):

    prank0 = []
    prank1 = []
    prank2 = []
    prank3 = []
    prank4 = []
    prank5 = []

    Rp = Zp(p)
    vp = Rp.valuation()
    S.<y> = Zp(p, prec = 100, type = 'capped-rel', print_mode = 'series')[]

    for polyns in l:
        if polyns.is_irreducible():
            d = polyns.degree()
            vpa5 = vp(polyns[d - 5])
            vpa4 = vp(polyns[d - 4])
            vpa3 = vp(polyns[d - 3])
            vpa2 = vp(polyns[d - 2])
            vpa1 = vp(polyns[d - 1])
        
            if vpa5 == 0:
                prank5.append(polyns)
            else:
                rootsof = polyns.roots(Rp)
                rv = []
                if not len(rootsof) == 0:
                    for roots in rootsof:
                        rv.append(vp(roots[0]))
            
                coeffsinZp = []
                for co in range(polyns.degree() + 1):
                    coeffsinZp.append(polyns[co])
                    polynsinZp = S(coeffsinZp)
            
                factorinZp = polynsinZp.factor()
                fd = []
                for fa in factorinZp:
                    fd.append(fa[0].degree())
            
                if vpa5 >= n/2 and vpa4 == 0 and n/2 not in rv:
                    prank4.append(polyns)
            
                elif vpa5 >= n and vpa4 >= n/2 and vpa3 == 0 and n/2 not in rv:
                    prank3.append(polyns)
            
                elif vpa5 == n and vpa4 >= n*2/3 and vpa3 >= n/3 and vpa2 == 0 and (n*2/3 not in rv) and (n/3 not in rv):
                    prank2.append(polyns)
            
                elif vpa5 >= n*3/2 and vpa4 >= n and vpa3 >= n/2 and vpa2 == 0 and (n/2 not in rv) and 3 not in fd:
                    prank2.append(polyns)            

                elif vpa5 == n and vpa4 >= 3*n/4 and vpa3 >= n/2 and vpa2 >= n/4 and vpa1 == 0 and (n/4 not in rv) and (n*3/4 not in rv) and 2 not in fd:
                    prank1.append(polyns) 

                elif vpa5 >= n*3/2 and vpa4 == n and vpa3 >= 2*n/3 and vpa2 >= n/3 and vpa1 == 0 and (n/3 not in rv) and (n*2/3 not in rv) and (n/2 not in rv):
                    prank1.append(polyns) 
                    
                elif vpa5 >= 2*n and vpa4 >= n*3/2 and vpa3 >= n and vpa2 >= n/2 and vpa1 == 0 and (n/2 not in rv) and 3 not in fd:
                    prank1.append(polyns) 
                    
                elif vpa5 == n and vpa4 >= n*4/5 and vpa3 >= n*3/5 and vpa2 >= n*2/5 and vpa1 >= n/5 and len(rv) == 0 and 2 not in fd:
                    prank0.append(polyns) 
                
                elif vpa5 >= n*3/2 and vpa4 == n and vpa3 >= n*3/4 and vpa2 >= n/2 and vpa1 >= n/4 and len(rv) == 0 and fd.count(2) == 1:
                    prank0.append(polyns)
                
                elif vpa5 >= 2*n and vpa4 >= n*3/2 and vpa3 == n and vpa2 >= n*2/3 and vpa1 >= n/3 and len(rv) == 0:
                    prank0.append(polyns)

                elif vpa5 == 2*n and vpa4 >= n*8/5 and vpa3 >= n*6/5 and vpa2 >= n*4/5 and vpa1 >= n*2/5 and len(rv) == 0 and (2 not in fd):
                    prank0.append(polyns)
                    
                elif vpa5 >= n*5/2 and vpa4 >= 2*n and vpa3 >= n*3/2 and vpa2 >= n and vpa1 >= n/2 and len(rv) == 0 and (3 not in fd) and (5 not in fd):
                    prank0.append(polyns)
                
        elif len(polyns.factor()) == 1:
            if n % 5 == 0 and polyns.factor()[0][1] == 5:
                omega = polyns.factor()[0][0][1]
                if abs(omega) < 2 * sqrt(q) and omega % q^(1/5) == 0:
                    if omega % q^(2/5) == 0:
                        if not omega / q^(2/5) % p == 0:
                            prank0.append(polyns)
                    elif not omega / q^(1/5) % p == 0:
                        prank0.append(polyns)

    total = prank0 + prank1 + prank2 + prank3 + prank4 + prank5
    
    return prank5, prank4, prank3, prank2, prank1, prank0, total

prank5, prank4, prank3, prank2, prank1, prank0, total = CharpolysofsimpleAVofdim5frombounds(l,p,n)

#beneath is the comparison of characteristic polynomials of simple abelian varieties of dimension 5 over F_q with data from LMFDB


# In[132]:


Missing5 = []

for m in charpolydim5prank5:
    if m not in prank5:
        Missing5.append(m)

extra5 = []

for M in prank5:
    if M not in charpolydim5prank5:
        extra5.append(M)

print("number of characteristic polynomials p-rank 5 in LMFDB:", len(data_dim5_prank5))
print("number of characteristic polynomials p-rank 5 according to code:", len(prank5))
print("characteristic polynomials in LMFDB but not in code:", Missing5)
print("characteristic polynomials in code, but not in LMFDB:", extra5)


# In[126]:


Missing4 = []

for m in charpolydim5prank4:
    if m not in prank4:
        Missing4.append(m)

extra4 = []

for M in prank4:
    if M not in charpolydim5prank4:
        extra4.append(M)

print("number of characteristic polynomials p-rank 4 in LMFDB:", len(data_dim5_prank4))
print("number of characteristic polynomials p-rank 4 according to code:", len(prank4))
print("characteristic polynomials in LMFDB but not in code:", Missing4)
print("characteristic polynomials in code, but not in LMFDB:", extra4)


# In[127]:


Missing3 = []

for m in charpolydim5prank3:
    if m not in prank3:
        Missing3.append(m)

extra3 = []

for M in prank3:
    if M not in charpolydim5prank3:
        extra3.append(M)

print("number of characteristic polynomials p-rank 3 in LMFDB:", len(data_dim5_prank3))
print("number of characteristic polynomials p-rank 3 according to code:", len(prank3))
print("characteristic polynomials in LMFDB but not in code:", Missing3)
print("characteristic polynomials in code, but not in LMFDB:", extra3)


# In[128]:


Missing2 = []

for m in charpolydim5prank2:
    if m not in prank2:
        Missing2.append(m)

extra2 = []

for M in prank2:
    if M not in charpolydim5prank2:
        extra2.append(M)

print("number of characteristic polynomials p-rank 2 in LMFDB:", len(data_dim5_prank2))
print("number of characteristic polynomials p-rank 2 according to code:", len(prank2))
print("characteristic polynomials in LMFDB but not in code:", Missing2)
print("characteristic polynomials in code, but not in LMFDB:", extra2)


# In[129]:


Missing1 = []

for m in charpolydim5prank1:
    if m not in prank1:
        Missing1.append(m)

extra1 = []

for M in prank1:
    if M not in charpolydim5prank1:
        extra1.append(M)

print("number of characteristic polynomials p-rank 1 in LMFDB:", len(data_dim5_prank1))
print("number of characteristic polynomials p-rank 1 according to code:", len(prank1))
print("characteristic polynomials in LMFDB but not in code:", Missing1)
print("characteristic polynomials in code, but not in LMFDB:", extra1)


# In[130]:


Missing0 = []

for m in charpolydim5prank0:
    if m not in prank0:
        Missing0.append(m)

extra0 = []

for M in prank0:
    if M not in charpolydim5prank0:
        extra0.append(M)

print("number of characteristic polynomials p-rank 0 in LMFDB:", len(data_dim5_prank0))
print("number of characteristic polynomials p-rank 0 according to code:", len(prank0))
print("characteristic polynomials in LMFDB but not in code:", Missing0)
print("characteristic polynomials in code, but not in LMFDB:", extra0)


# In[133]:


Missingtotal = Missing0 + Missing1 + Missing2 + Missing3 + Missing4 + Missing5
extratotal = extra0 + extra1 + extra2 + extra3 + extra4 + extra5

print("total number of characteristic polynomials in LMFDB:", len(charpolydim5))
print("total number of characteristic polynomial according to code:", len(total))
print("characteristic polynomials in LMFDB but not in code:", Missingtotal)
print("characteristic polynomials in code but not in LMFDB:", extratotal)


# In[ ]:




