def find_line_data(ionline):
    
    ## this needs to be more sophisticated 
    ## Morton 2003 
    ionline = ionline.replace(' ', '')
    line_struct = {'SiIV1393': {'wave': 1393.7602, 'fval': 5.13E-01, 'Ref': 'Morton+2003'}, 
                   'SiIV1402': {'wave': 1402.7729, 'fval': 2.54E-01, 'Ref': 'Morton+2003'}, 
                   'CIV1548': {'wave': 1548.204, 'fval': 1.899E-01, 'Ref': 'Morton+2003'}, 
                   'CIV1550': {'wave': 1550.781, 'fval': 9.475E-02, 'Ref': 'Morton+2003'}, 
                   'SII1250': {'wave': 1250.578, 'fval': 5.43E-03, 'Ref': 'Morton+2003'}, 
                   'SII1253': {'wave': 1253.805, 'fval': 1.09E-02, 'Ref': 'Morton+2003'}, 
                   'SII1259': {'wave': 1259.518, 'fval': 1.66E-02, 'Ref': 'Morton+2003'}
                  }
    return line_struct[ionline]
