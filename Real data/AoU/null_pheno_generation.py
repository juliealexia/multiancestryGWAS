for eth in ['AFR','AMR','EAS','EUR','MID','SAS']:
    file_name = f'{eth}_pca.eigenvec'
    pcs_local = pd.read_csv(file_name, sep='\t')
    file_name2 = f'{eth}_pca.eigenval'
    eigenval_local = pd.read_csv(file_name2, sep='\t',header=None)[0]
    pcs_local['IID'] = pcs_local['IID'].astype(int)
    pheno = pd.DataFrame({'FID': [0]*pcs_local.shape[0],'IID': pcs_local['IID']})
    #Set % of variance captured by local PCs 
    for tarvarlocal in [0,0.002,0.005,0.01,0.02,0.05]:
        new_col_name = f"top10_Null_{int(tarvarlocal*1000)}L"
        varpc1 = pcs_local['PC1'].var()
        w1 = eigenval_local[0]/eigenval_local.sum()
        beta1 = np.sqrt((tarvarlocal*w1)/varpc1)
        varpc2 = pcs_local['PC2'].var()
        w2 = eigenval_local[1]/eigenval_local.sum()
        beta2 = np.sqrt((tarvarlocal*w2)/varpc2)
        varpc3 = pcs_local['PC3'].var()
        w3 = eigenval_local[2]/eigenval_local.sum()
        beta3 = np.sqrt((tarvarlocal*w3)/varpc3)
        varpc4 = pcs_local['PC4'].var()
        w4 = eigenval_local[3]/eigenval_local.sum()
        beta4 = np.sqrt((tarvarlocal*w4)/varpc4)
        varpc5 = pcs_local['PC5'].var()
        w5 = eigenval_local[4]/eigenval_local.sum()
        beta5 = np.sqrt((tarvarlocal*w5)/varpc5)
        varpc6 = pcs_local['PC6'].var()
        w6 = eigenval_local[5]/eigenval_local.sum()
        beta6 = np.sqrt((tarvarlocal*w6)/varpc6)
        varpc7 = pcs_local['PC7'].var()
        w7 = eigenval_local[6]/eigenval_local.sum()
        beta7 = np.sqrt((tarvarlocal*w7)/varpc7)
        varpc8 = pcs_local['PC8'].var()
        w8 = eigenval_local[7]/eigenval_local.sum()
        beta8 = np.sqrt((tarvarlocal*w8)/varpc8)
        varpc9 = pcs_local['PC9'].var()
        w9 = eigenval_local[8]/eigenval_local.sum()
        beta9 = np.sqrt((tarvarlocal*w9)/varpc9)
        varpc10 = pcs_local['PC10'].var()
        w10 = eigenval_local[9]/eigenval_local.sum()
        beta10 = np.sqrt((tarvarlocal*w10)/varpc10)
        noise = np.random.normal(0, np.sqrt(1-tarvarlocal), pcs_local.shape[0])
        pheno[new_col_name] = beta1 * pcs_local['PC1'] + beta2 * pcs_local['PC2'] + beta3 * pcs_local['PC3'] + beta4 * pcs_local['PC4'] + beta5 * pcs_local['PC5'] + beta6 * pcs_local['PC6'] + beta7 * pcs_local['PC7'] + beta8 * pcs_local['PC8'] + beta9 * pcs_local['PC9'] + beta10 * pcs_local['PC10'] + noise

    out_name = f'{eth}_NULL_first10.txt'
    pheno.to_csv(out_name, sep='\t', index=False)


pcs_global = pd.read_csv('ALL_pca.eigenvec', sep='\t')
#Rename PC1-PC10 as PC1_G-PC10_G to indicate global PCs
cols = pcs_global.columns.tolist()
new_cols = cols[:2] + [col + "_G" for col in cols[2:12]] + cols[12:]
pcs_global.columns = new_cols
#Remove FID column as all zeros and easier for merging later on
pcs_global = pcs_global.drop(pcs_global.columns[0], axis=1)
eigenval_global = pd.read_csv('ALL_pca.eigenval', sep='\t',header=None)[0]
np.random.seed(3007)
pheno = pd.DataFrame({'FID': [0]*pcs_global.shape[0],'IID': pcs_global['IID']})
    #Set % of variance captured by global PCs 
    for tarvarlocal in [0,0.002,0.005,0.01,0.02,0.05]:
        new_col_name = f"top10_Null_{int(tarvarlocal*1000)}L"
        varpc1 = pcs_global['PC1'].var()
        w1 = eigenval_global[0]/eigenval_global.sum()
        beta1 = np.sqrt((tarvarlocal*w1)/varpc1)
        varpc2 = pcs_global['PC2'].var()
        w2 = eigenval_global[1]/eigenval_global.sum()
        beta2 = np.sqrt((tarvarlocal*w2)/varpc2)
        varpc3 = pcs_global['PC3'].var()
        w3 = eigenval_global[2]/eigenval_global.sum()
        beta3 = np.sqrt((tarvarlocal*w3)/varpc3)
        varpc4 = pcs_global['PC4'].var()
        w4 = eigenval_global[3]/eigenval_global.sum()
        beta4 = np.sqrt((tarvarlocal*w4)/varpc4)
        varpc5 = pcs_global['PC5'].var()
        w5 = eigenval_global[4]/eigenval_global.sum()
        beta5 = np.sqrt((tarvarlocal*w5)/varpc5)
        varpc6 = pcs_global['PC6'].var()
        w6 = eigenval_global[5]/eigenval_global.sum()
        beta6 = np.sqrt((tarvarlocal*w6)/varpc6)
        varpc7 = pcs_global['PC7'].var()
        w7 = eigenval_global[6]/eigenval_global.sum()
        beta7 = np.sqrt((tarvarlocal*w7)/varpc7)
        varpc8 = pcs_global['PC8'].var()
        w8 = eigenval_global[7]/eigenval_global.sum()
        beta8 = np.sqrt((tarvarlocal*w8)/varpc8)
        varpc9 = pcs_global['PC9'].var()
        w9 = eigenval_global[8]/eigenval_global.sum()
        beta9 = np.sqrt((tarvarlocal*w9)/varpc9)
        varpc10 = pcs_global['PC10'].var()
        w10 = eigenval_global[9]/eigenval_global.sum()
        beta10 = np.sqrt((tarvarlocal*w10)/varpc10)
        noise = np.random.normal(0, np.sqrt(1-tarvarlocal), pcs_global.shape[0])
        pheno[new_col_name] = beta1 * pcs_global['PC1'] + beta2 * pcs_global['PC2'] + beta3 * pcs_global['PC3'] + beta4 * pcs_global['PC4'] + beta5 * pcs_global['PC5'] + beta6 * pcs_global['PC6'] + beta7 * pcs_global['PC7'] + beta8 * pcs_global['PC8'] + beta9 * pcs_global['PC9'] + beta10 * pcs_global['PC10'] + noise

    out_name = 'ALL_NULL_first10G.txt'
    pheno.to_csv(out_name, sep='\t', index=False)

#Subset for individual files by ancestry