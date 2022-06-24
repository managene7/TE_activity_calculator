def RT_similarity(csv_name):
    import csv
    input_csv=csv.reader(open(csv_name,'r'))
    pair_list=[]
    idty_dic={}
    
    for row in input_csv:
        if row[0] != 'Query_name':
            #if '_continued' not in row[4]:
            if row[0].split('_continued')[0] != row[4].split('_continued')[0]:
                pair_list.append([float(row[10]),row[0].split('_continued')[0],row[4].split('_continued')[0]])
                if row[0].split('_continued')[0] not in idty_dic:
                    idty_dic[row[0].split('_continued')[0]]={row[4].split('_continued')[0]:float(row[10])}
                else:
                    idty_dic[row[0].split('_continued')[0]][row[4].split('_continued')[0]]=float(row[10])
    #print idty_dic
    pair_list.sort()
    pair_list.reverse()
    
    reference=[]
    identity_list=[]
    zero_div=0
    phylo_list=[]
    for cont in pair_list:
        init=0
        for check in reference:
            if cont[1] in check and cont[2] in check:
                init=1
        if init==0:
            left_list=[]
            right_list=[]
            
            n=0
            for check in reference:        
                if cont[1] in check:
                    left_list.append(check)
                    l_posi=n
                n=n+1     
    
            k=0    
            for check in reference:
                if cont[2] in check:
                    right_list.append(check)
                    r_posi=k
                k=k+1
            if len(left_list) >=2:
                print("\n\nError occurred in left list..\n\n")
                quit()
            if len(right_list) >=2:
                print("\n\nError occurred in right list..\n\n")
                quit()
            if len(left_list)==0 and len(right_list)==0:
                reference.append([cont[1],cont[2]])
                
                try:
                    identity_list.append(idty_dic[cont[1]][cont[2]])
                    #print [idty_dic[cont[1]][cont[2]],[cont[1],cont[2]]]
                    phylo_list.append([idty_dic[cont[1]][cont[2]],[cont[1],cont[2]]])###############
                    
                except:
                    print('There are not aligned sequences..00')
                    
                
            elif len(left_list)==1 and len(right_list)==0:
                tot_idty=0
                count=0
                temp_list=[]
                #print left_list[0]
                for element in left_list[0]:
                    if element in idty_dic:
                        if cont[2] in idty_dic[element]:
                            temp_list.append(element)
                            count=count+1
                            tot_idty=tot_idty+idty_dic[element][cont[2]]
                        
                    else:
                        print('There are not aligned sequences..10')
                try:
                    #print count
                    avg=tot_idty/float(count)
                    reference[l_posi].append(cont[2])
                    identity_list.append(avg)
                    
                    phylo_list.append([avg,[temp_list,cont[2]]])#####################
                except:
                    zero_div=zero_div+1
                    print("There was 0 division..%s" % str(zero_div))
                    pass
                
            elif len(left_list)==0 and len(right_list)==1:
                tot_idty=0
                count=0
                temp_list=[]
                for element in right_list[0]:
                    if element in idty_dic:
                        if cont[1] in idty_dic[element]:
                            temp_list.append(element)
                            count=count+1
                            tot_idty=tot_idty+idty_dic[element][cont[1]]
                    else:
                        print('There are not aligned sequences..01')
                try:
                    #print count
                    avg=tot_idty/float(count)
                    
                    reference[r_posi].append(cont[1])
                    identity_list.append(avg)
                    
                    phylo_list.append([avg,[cont[1],temp_list]])#####################
                except:
                    zero_div=zero_div+1
                    print("There was 0 division..%s" % str(zero_div))
                    pass
                
            elif len(left_list)==1 and len(right_list)==1:
                tot_idty=0
                count=0
                left_data=left_list[0]
                right_data=right_list[0]
                temp_list1=[]
                temp_list2=[]
                for element1 in left_data:
                    for element2 in right_data:
                        if element1 in idty_dic:
                            temp_list1.append(element1)
                            if element2 in idty_dic[element1]:
                                temp_list2.append(element2)
                                #print element1, element2
                                count=count+1
                                tot_idty=tot_idty+idty_dic[element1][element2]
                            else:
                                print('There are not aligned sequences..11')
                
                try:
                    #print count
                    avg=tot_idty/float(count)
                    
                    identity_list.append(avg)
                    
                    phylo_list.append([avg,[list(set(temp_list1)), list(set(temp_list2))]])#####################
                except:
                    zero_div=zero_div+1
                    print("There was 0 division..%s" % str(zero_div))
                    pass
                new_list=left_list[0]+right_list[0]
                reference.remove(left_list[0])
                reference.remove(right_list[0])
                
                reference.append(new_list)
    identity_list.sort()
    identity_list.reverse()
    newic_list=[]
    for cont in phylo_list:
        #print (cont)
        #cont[1].sort()
        name_list=[]
        newic_sub=''
        for sp in cont[1]:
            #print sp
            if type(sp)!=list:
                name_list.append(sp)
                #print sp
                newic_sub=newic_sub+sp+":"+str(1.0-float(cont[0]))+','
            else:
                sub_name_list=[]
                for sp2 in sp:
                    name_list.append(sp2)
                    sub_name_list.append(sp2)
                sub_name_list.sort()
                for newic_cont in newic_list:
                    if sub_name_list in newic_cont:
                        #newic_sub=newic_sub+newic_cont[1]+str(float(newic_cont[2])-float(cont[0]))+':'+str(cont[0])+','
                        newic_sub=newic_sub+newic_cont[1]+':'+str((1.0-float(cont[0]))-newic_cont[2])+','
        newic_sub=newic_sub[:-1]    
        name_list.sort()
        #print name_list,"("+newic_sub+")"
        newic_list.append([name_list,"("+newic_sub+")",1-float(cont[0])])
        
    out_newic=open(csv_name[:-4]+'_newic_for_phylo.nwk','w')
    out_newic.write(newic_list[-1][1]+';')
    #print newic_list[-1][1]+';'
        #print cont
    return identity_list


#!/usr/bin/env python

#___this code needs '0_blast_xml_parsing_memory_safe_v2.4_sequence_added_in_csv_position_marked_split_included.py' file_

import csv

#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-out':'RT_activity'}
for i in range(len(args)):
    if args[i].startswith("-"):
        try:
            option_dict[args[i]]=args[i+1]
        except:
            if args[0]=="-help":
                print("""
_____________________________________________________________________________

Usage;

-csv        csv file name of blast result (default is 'RT_activity')
_____________________________________________________________________________
""")
                quit()

csv_out=csv.writer(open(option_dict['-csv']+"_UPGMA_distance.csv",'w', newline=""))
value_list=RT_similarity(option_dict['-csv'])
cl_name=[option_dict['-csv']]
cont=cl_name+value_list
csv_out.writerow(cont)

            
            
            
            
            
            

