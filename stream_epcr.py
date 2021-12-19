#!/usr/bin/env python
# coding: utf-8

# In[ ]:

from __future__ import division
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from urllib.error import URLError
import math
import requests
import io
pd.set_option('precision', 1)

st.set_page_config(layout="wide")

url = "https://raw.githubusercontent.com/j0n0curry/ePCR_viewer/master/ALT8_concord1.csv"
download = requests.get(url).content

# Reading the downloaded content and turning it into a pandas dataframe

concordance = pd.read_csv(url)



#percentiles
def Q25(x):
    return x.quantile(0.25)

def Q50(x):
    return x.quantile(0.5)

# 90th Percentile
def Q75(x):
    return x.quantile(0.75)



stats_FAM = concordance.groupby(['Sample','CALL_N3'])['FAM_N3'].agg(['count', 'mean','std','min',Q25, Q50, Q75, 'max'])
#print(stats_FAM)
print('-'*30)

CI95_hi_FAM = []
CI95_lo_FAM = []
CV_run_FAM = []

for i in stats_FAM.index:
    c,m,s,t,u,q1,q2,v = round(stats_FAM.loc[i])
    CI95_hi_FAM.append(m + 1.95*s/math.sqrt(c))
    CI95_lo_FAM.append(m - 1.95*s/math.sqrt(c))
    #CV_run_FAM.append(s/m*100)

stats_FAM['ci95_lo_FAM'] = CI95_lo_FAM
stats_FAM['ci95_hi_FAM'] = CI95_hi_FAM
#stats_FAM['CV%_FAM'] = CV_run_FAM

print(stats_FAM)


print('-'*30)

print('nomralised_FAM_N1N2 by sample type')

stats_nFAM = concordance.groupby(['Sample','CALL_N3'])['n_FAM'].agg(['count', 'mean','min', 'std',Q25, Q50, Q75, 'max'])

print('-'*30)

CI95_hi_nFAM = []
CI95_lo_nFAM = []
CV_run_nFAM = []

for i in stats_nFAM.index:
    c,m,s,t,u,q1,q2,v =(stats_nFAM.loc[i])
    CI95_hi_nFAM.append(m + 1.95*s/math.sqrt(c))
    CI95_lo_nFAM.append(m - 1.95*s/math.sqrt(c))
    CV_run_nFAM.append(s/m*100)

stats_nFAM['Confidence Interval 95% low nFAM'] = CI95_lo_nFAM
stats_nFAM['ci95_hi_nFAM'] = CI95_hi_nFAM
stats_nFAM['CV%_nFAM'] = CV_run_nFAM
#stats_nFAM['%Percent_detected'] = result['N1N2_detected'] / TOT*100
print(stats_nFAM)


fig2b = px.scatter(concordance, x= 'n_VIC_N3', y = 'n_FAM_N3',color = 'Sample')

#fig2b.show()
#fig2b.write_html('Concordance_N3.html')


fig1bbnbb = px.scatter(concordance, x= 'Order', y = 'n_VIC_N3', color = 'CONCORD')
fig1bbnbb.update_yaxes(range=[0, 6])
#fig1bbnbb.show()
#fig1bbnbb.write_html('concordance_N3__monitor_normRNaseP.html')

figROX = px.scatter(concordance, x= 'Order', y = 'ROX_N3', color = 'CALL_N3', title = 'Nexar 3 Dispense Trace ROX')
figROX.update_yaxes(gridwidth = 0.0002, gridcolor ='red')

figROX.add_trace(go.Scatter(
    x=[concordance.Order.min(), concordance.Order.min()],
    y=[1600, 1600],
    mode="lines",
    name="1600  RFU Lower Cutoff Limit",
    text=["LCL"],
    #text=["ROX 1600 lower cutoff"],
    textposition="top center",
    line=dict(color="grey")
))
#figrp.show()


figrp = px.scatter(concordance, x= 'ROX_N3', y = 'FAM_N3' ,color = 'CONCORD')
figrp.update_xaxes(range=[1000, 6000])
figrp.update_yaxes(range=[0, 50000])


figrp.add_trace(go.Scatter(
    x=[1600, 1600],
    y=[50000, -100],
    mode="lines",
    name="1600  RFU Lower Cutoff Limit",
    text=["LCL"],
    #text=["ROX 1600 lower cutoff"],
    textposition="top center",
    line=dict(color="grey")
))
#figrp.show()

print('CFO_RFU by sample type')

stats_CFO = concordance.groupby(['Sample','CALL_N3'])['VIC_N3'].agg(['count', 'mean','std','min',Q25, Q50, Q75, 'max'])

CI95_hi_CFO = []
CI95_lo_CFO = []
CV_run_CFO = []

for i in stats_CFO.index:
    c,m,s,t,u,q1,q2,v = round(stats_CFO.loc[i])
    CI95_hi_CFO.append(m + 1.95*s/math.sqrt(c))
    CI95_lo_CFO.append(m - 1.95*s/math.sqrt(c))
    CV_run_CFO.append(s/m*100)

stats_CFO['ci95_lo_CFO'] = CI95_lo_CFO
stats_CFO['ci95_hi_CFO'] = CI95_hi_CFO
stats_CFO['CV%_CFO'] = CV_run_CFO

print(stats_CFO)




print('-'*30)

print('normalised_CFO_RNaseP by sample type')

stats_nCFO = concordance.groupby(['Sample','CALL_N3'])['n_VIC_N3'].agg(['count', 'mean','std','min',Q25, Q50, Q75, 'max'])
print(stats_nCFO)
print('-'*30)

CI95_hi_nCFO = []
CI95_lo_nCFO = []
CV_run_nCFO = []

for i in stats_nCFO.index:
    c,m,s,t,u,q1,q2,v = (stats_nCFO.loc[i])
    CI95_hi_nCFO.append(m + 1.95*s/math.sqrt(c))
    CI95_lo_nCFO.append(m - 1.95*s/math.sqrt(c))
    CV_run_nCFO.append(s/m*100)

stats_nCFO['ci95_lo_nCFO'] = CI95_lo_nCFO
stats_nCFO['ci95_hi_nCFO'] = CI95_hi_nCFO
stats_nCFO['CV%_nCFO'] = CV_run_nCFO

print(stats_nCFO)


figN1 = px.scatter(concordance, x= 'Order', y = 'n_FAM_N3' ,color = 'CALL_N3', title = 'Nexar 3 - N1 N2 Calls')

figN1.add_trace(go.Scatter(
    y=[10, 10],
    x=[concordance.Order.min(), concordance.Order.max()],
    mode="lines+markers+text",
    name="Lower_10_Positive_Boundary",
    text=["10"],
    #text=["ROX 1600 lower cutoff"],
    textposition="top center",
    line=dict(color="red")
))



figN1.add_trace(go.Scatter(
     y=[9, 9],
     x=[concordance.Order.min(), concordance.Order.max()],
     mode="lines+markers+text",
     name="Lower_9_Positive_Boundary",
     text=["9"],
     textposition="top center"))



figN1.update_xaxes(showgrid = True, gridwidth = 0.0002, gridcolor = 'grey')
figN1.update_yaxes(range=[0, 16],showgrid = True, gridwidth = 0.0002, gridcolor = 'grey')


#figN1.show()


figcomp = px.scatter(concordance, x= 'FAM', y = 'FAM_N3' ,color = 'CONCORD', title = 'Concordance')
figcomp.update_xaxes(range=[0, 50000], showgrid = True, gridwidth = 0.0002, gridcolor = 'grey')
figcomp.update_yaxes(range=[0, 50000], showgrid = True, gridwidth = 0.0002, gridcolor = 'grey')


figcomp.add_trace(go.Scatter(
    x=[0, 50000],
    y=[0, 50000],
    mode="lines",
    name="Concordance line",
    text=["Concordance line"],
    #text=["Concordance line"],
    textposition="top center",
    line=dict(color="red")
))
#figrp.show()

# Plot!
st.title('ePCR viewer')

col1, col2 = st.columns(2)

with col1:
	st.plotly_chart(figrp, use_container_width=True)
with col2:
	st.plotly_chart(fig2b, use_container_width=True)

st.plotly_chart(figROX,use_container_width = True)
st.plotly_chart(figN1, use_container_width=True)
st.plotly_chart(fig1bbnbb, use_container_width=True)
st.plotly_chart(figcomp, use_container_width = True)


st.table(stats_nFAM)
# Streamlit widgets automatically run the script from top to bottom. Since
# this button is not connected to any other logic, it just causes a plain
# rerun.
st.button("Re-run")