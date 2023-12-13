# -*- coding: utf-8 -*-
"""
Extract data from GSEA analysis output

"""

from bs4 import BeautifulSoup
import requests
import pandas as pd
import os

gsea_dir=os.listdir("C:\\Users\\jxf43\\gsea_home\\output\\jul23")
gsea_dir_name="C:\\Users\\jxf43\\gsea_home\\output\\jul23"

for dir_name in gsea_dir:
    dir_path_name=os.path.join(gsea_dir_name,dir_name)
    files=os.listdir(dir_path_name)
    result=[]
    for f in files:
        if f.startswith("GOCC_") and f.endswith("html"):
            filename=os.path.join(dir_path_name,f)
            with open(filename, 'r') as file:
                html_content = file.read()
            table_data=[]  
            soup=BeautifulSoup(html_content, 'html.parser')
            table = soup.find('table')
            if table:
                # Find all rows (<tr>) in the table.
                rows = table.find_all('tr')

                for row in rows:
                    # Find all cells (<td> or <th>) in each row.
                    cells = row.find_all(['td', 'th'])
                    
                    # Extract the text from each cell and append it to the table_data list.
                    row_data = [cell.get_text(strip=True) for cell in cells]
                    table_data.append(row_data)
                df = pd.DataFrame(table_data)
                df=df.T
                df.columns = df.iloc[0]
                df = df.iloc[1:]
                #geneset=str(df['GeneSet'])
                result.append(df)
    combined_df = pd.concat(result, ignore_index=True)
    file_name=dir_name+"_gsea_result.csv"
    combined_df.to_csv(file_name,index=False)


dir_path="C:\\Users\\jxf43\\gsea_home\\output\\jul21\\my_analysis.Gsea.1689972025057"
files=os.listdir(dir_path)
result=[]
for f in files:
    if f.startswith("GOCC_") and f.endswith("html"):
        filename=os.path.join(dir_path,f)
        with open(filename, 'r') as file:
            html_content = file.read()
        table_data=[]  
        soup=BeautifulSoup(html_content, 'html.parser')
        table = soup.find('table')
        if table:
            # Find all rows (<tr>) in the table.
            rows = table.find_all('tr')

            for row in rows:
                # Find all cells (<td> or <th>) in each row.
                cells = row.find_all(['td', 'th'])
                
                # Extract the text from each cell and append it to the table_data list.
                row_data = [cell.get_text(strip=True) for cell in cells]
                table_data.append(row_data)
            df = pd.DataFrame(table_data)
            df=df.T
            df.columns = df.iloc[0]
            df = df.iloc[1:]
            #geneset=str(df['GeneSet'])
            result.append(df)

combined_df = pd.concat(result, ignore_index=True)
combined_df.to_csv("local_combined_df.csv",index=False)


