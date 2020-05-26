#
#   This is a module file for b_samples.py
#   You need to define the dictionaries of transition data
#   For now, the key is case name (Mg26) and the entry is a list of dictionaries
#   Fox 2019

def get_trans_data(case_name,op_type):
    # transition data. trans[<case>] is a list of dictionaries with keys 2Ji, ni, 2Jf, nf , B
    # 2Ji and 2Jf are initial and final 2J
    # ni and nf are 1,2,3,... corresponding to the first, second, third, ... state with that J-value

    trans = {}
    trans['Mg26'] =[
            {'type':'E2','2Ji': 4  ,'ni': 1   ,'2Jf': 0,'nf': 1 ,'B':0.0},
            {'type':'E2','2Ji': 4  ,'ni': 2   ,'2Jf': 0,'nf': 1 ,'B':0.0},
            {'type':'E2','2Ji': 0  ,'ni': 2   ,'2Jf': 4,'nf': 1 ,'B':0.0},
            {'type':'E2','2Ji': 4  ,'ni': 4   ,'2Jf': 0,'nf': 1 ,'B':0.0}]
    trans['F18'] =[
            {'type':'M1','2Ji': 0, 'ni': 1, '2Jf': 2, 'nf': 1, 'B': 0.0},
            {'type':'M1','2Ji': 2, 'ni': 2, '2Jf': 2, 'nf': 1, 'B': 0.0},
            {'type':'M1','2Ji': 2, 'ni': 2, '2Jf': 0, 'nf': 1, 'B': 0.0},
            {'type':'M1','2Ji': 6, 'ni': 2, '2Jf': 6, 'nf': 1, 'B': 0.0},
            {'type':'M1','2Ji': 6, 'ni': 2, '2Jf': 4, 'nf': 1, 'B': 0.0},
            {'type':'M1','2Ji': 4, 'ni': 3, '2Jf': 4, 'nf': 1, 'B': 0.0},
            {'type':'M1','2Ji': 4, 'ni': 3, '2Jf': 6, 'nf': 1, 'B': 0.0},
            {'type':'M1','2Ji': 4, 'ni': 3, '2Jf': 2, 'nf': 2, 'B': 0.0},
            {'type':'M1','2Ji': 6, 'ni': 3, '2Jf': 4, 'nf': 2, 'B': 0.0},
            {'type':'M1','2Ji': 2, 'ni': 4, '2Jf': 4, 'nf': 2, 'B': 0.0},
            {'type':'M1','2Ji': 8, 'ni': 2, '2Jf': 8, 'nf': 1, 'B': 0.0}]
    trans['Al26'] =[
            {'type':'E2','2Ji': 6  ,'ni': 1   ,'2Jf': 10,'nf': 1 ,'B':0.0},
            {'type':'E2','2Ji': 2  ,'ni': 2   ,'2Jf': 6,'nf': 1 ,'B':0.0},
            {'type':'E2','2Ji': 4  ,'ni': 2   ,'2Jf': 0,'nf': 1 ,'B':0.0},
            {'type':'E2','2Ji': 2  ,'ni': 3   ,'2Jf': 6,'nf': 1 ,'B':0.0},
            {'type':'E2','2Ji': 6  ,'ni': 2   ,'2Jf': 10,'nf': 1 ,'B':0.0},
            {'type':'E2','2Ji': 6  ,'ni': 2   ,'2Jf': 2,'nf': 1 ,'B':0.0},
            {'type':'M1','2Ji': 2  ,'ni': 1   ,'2Jf': 0,'nf': 1 ,'B':0.0},
            {'type':'M1','2Ji': 2  ,'ni': 2   ,'2Jf': 0,'nf': 1 ,'B':0.0},
            {'type':'M1','2Ji': 2  ,'ni': 3   ,'2Jf': 0,'nf': 1 ,'B':0.0},
            {'type':'M1','2Ji': 2  ,'ni': 4   ,'2Jf': 0,'nf': 1 ,'B':0.0},
            {'type':'M1','2Ji': 4  ,'ni': 5   ,'2Jf': 2,'nf': 1 ,'B':0.0},
            {'type':'M1','2Ji': 4  ,'ni': 5   ,'2Jf': 6,'nf': 1 ,'B':0.0},
            {'type':'M1','2Ji': 2  ,'ni': 5   ,'2Jf': 0,'nf': 1 ,'B':0.0}]
    trans['Si32'] =[
            {'type':'GT','2Ji': 0  ,'ni': 1   ,'2Jf': 2,'nf': 1 ,'B':0.0}]

    trans_pick = {}
#    trans_pick[case_name] = trans[case_name]
    trans_pick[case_name] = []
    for x in trans[case_name]:
        if x['type']==op_type:
            trans_pick[case_name].append(x)

    return trans_pick


