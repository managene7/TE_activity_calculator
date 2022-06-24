#Parsing of XML file that contains BLAST results


#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-lim_evalue':10, '-ov_len':30, '-seq_parse':2, '-lim_identity':0.6,'-hsp_identity':0.6,'-lim_score':150,'-hit_function':'1','-num_hsp':10000000,'-num_hit':10000000,'-out':""}
for i in range(len(args)):
    if args[i].startswith("-"):
        try:
            option_dict[args[i]]=args[i+1]
        except:
            if args[0]=="-help":
                print("""
________________________________________________________________________________


Usage;

If default value exists, corresponding value can be omitted.

-xml            file name of blast result (xml format)
-out            output file name without extension that you want to (ex.. test)
-seq_parse      when you want to write the matched sequence in csv file..
                1: write aligned sequence 2: don't write. Default is 2.
-lim_evalue     evalue threshold of BLAST hit. Default is e-10
-lim_identity   identity threshold of BLAST hit. Default is 0.60
-hsp_identity   identity threshold of BLAST hsp. Default is 0.60
-num_hit        number of hit to be parsed
-num_hsp        number of hsp to be parsed
-lim_score      score threshold of hsps in a hit. Default is 150.
-hit_function   to filter out 'unknown' or 'hypothetical' in hit description.
                1: not filtering 2: filtering. Default is 1
-ov_len         length threshold of result overlap result. Default is 30.

________________________________________________________________________________
""")
                quit()
#        else:
#            print "There maybe an option error!"
#            quit()






input_file=open(option_dict['-xml'], 'r')

#print "\n\nBLAST xml parsing in progression......" ##################################################################

import csv
import shelve
if option_dict['-out']=="":
    option_dict['-out']=option_dict['-xml'][:-4]+".csv"
csv_file=csv.writer(open("%s" % option_dict["-out"],"w", newline=""))
title_line=["Query_name","Query_length", "Query_from", "Query_to", "Hit_description", "Hit_length", "Hit_from", "Hit_to", "Score", "e-value", "Identity", "Query_coverage", "Subject_coverage"]
if int(option_dict['-seq_parse'])==1:
    title_line.append("query_seq")
    title_line.append("hit_seq")
csv_file.writerow(title_line)

#seq_name_list=open("%s.txt" % option_dict["-out"]+"-name_list", 'a') # Sequence name list file

query_def="temporary_query_def"
#This is MODULE 1 to read parameter of a blastp result by dictionary form
#structure {query_def:[query_len,[hit_num, hit_id, hit_len, total_score, E-value, mean_identity, query_coverage, subject_coverage]
#Pickle out_file will be returned


#________________________________________query circuit_____________________
parameter_dic={}
while 1:
    line=input_file.readline().strip()
    if not line:
        #print "\n\n\nBLAST xml parsing completed. Open the csv file by MicroSoft EXCEL..." ############################3
        break 

     
    if line.startswith("<Iteration>"):
        while 1:
            line=input_file.readline().strip()
            if line.startswith("<Iteration_query-def>"):
                query_def=line[21:-22]
            if line.startswith("<Iteration_query-len>"):
                query_len=line[21:-22]
                parameter_dic[query_def]=[int(query_len)]
                parameter_dic[query_def].append({})
            if line.startswith("<Iteration_hits></Iteration_hits>"):
                csv_file.writerow([query_def,parameter_dic[query_def][0], "-", "-", "-", "-", "-", "-", "-","-", "-", "-", "-", "-"]) 
                break
            if line.startswith("</Iteration>"):
                
                
#__________________garbage hit deletion______________________________



#This is MODULE 2 to analysis
    #Select maximum score hit
        #csv export

                if query_def != "temporary_query_def":


              #___________________________________hit deletion containing nameless___________________  

                    if str(option_dict['-hit_function'])=="2":
                        hit_id_list3=list(parameter_dic[query_def][1].keys())
        
                        hit_e_value_list=[]
                        for h_id in hit_id_list3:
                            hit_e_value_list.append([parameter_dic[query_def][1][h_id]['e-value'], (1000/parameter_dic[query_def][1][h_id]['total_score']), parameter_dic[query_def][1][h_id]['mean_identity'], h_id, parameter_dic[query_def][1][h_id]['hit_range_merged'], parameter_dic[query_def][1][h_id]['hit_def']])######### To delete Vitis vinifera hit, hit_def added
                        hit_e_value_list.sort()
        
        
                        hit_garbage_list=[]
                    
                        for key in hit_e_value_list:
                        
                            if 'hypothetical protein' in key[5]: ######### To delete nameless protein hit######
                                hit_garbage_list.append(key[3])######### To delete nameless protein hit######
                            elif 'predicted protein' in key[5]: ######### To delete nameless protein hit######
                                hit_garbage_list.append(key[3])######### To delete nameless protein hit######
                            elif 'uncharacterized' in key[5]:######### To delete nameless protein hit######
                                hit_garbage_list.append(key[3])######### To delete nameless protein hit######
                            elif 'uncharacterized protein' in key[5]:######### To delete nameless protein hit######
                                hit_garbage_list.append(key[3])######### To delete nameless protein hit######
                                
                            else:
                                break
           
                        if len(hit_garbage_list) >= 1:
                            for key in hit_garbage_list:
                                del parameter_dic[query_def][1][key]  #hit grabage deletion by e-value

                        hit_garbage_list=[] #To clear memory

            
                    #__________________garbage hit deletion by overlap______________________________
    
                    hit_id_list=list(parameter_dic[query_def][1].keys())
                    if len(hit_id_list) > 1:
                        hit_e_value_list=[]
                        for h_id in hit_id_list:
                            hit_e_value_list.append([                                
                                (1000/parameter_dic[query_def][1][h_id]['max_score']),
                                (100/parameter_dic[query_def][1][h_id]['identity_of_max_score']),
                                parameter_dic[query_def][1][h_id]['e-value'],
                                h_id,
                                parameter_dic[query_def][1][h_id]['hit_range_merged']])
                        hit_e_value_list.sort()

                        n=0
                        hit_garbage_list=[]
            
                        for key in hit_e_value_list:
                            n=n+1
                            if n>int(option_dict['-num_hit']):
                                hit_garbage_list.append(key[3])
                        
                            if n == len(hit_e_value_list):
                                break
                            for key_2 in hit_e_value_list[n:]:
                                if key[3] != key_2[3]:  # to avoid muplicated hit deletion
                                    if len(key_2[4].intersection(key[4])) > int(option_dict["-ov_len"]):
                                        hit_garbage_list.append(key_2[3])
                                else: pass


                        temp=set(hit_garbage_list)
                        hit_garbage_new=list(temp)

                        if len(hit_garbage_new) >= 1:
                            for key in hit_garbage_new:
                                del parameter_dic[query_def][1][key]  #hit grabage deletion by range comparison

                        hit_garbage_new=[] #To clear memory
                        temp=[] #To clear memory
                        hit_e_value_list=[] #To clear memory


              #___________________________________hit grabage deletion by e-value and identity threshold___________________  

                    hit_id_list2=list(parameter_dic[query_def][1].keys())
            
        
                    hit_e_value_list=[]
                    for h_id in hit_id_list2:
                        hit_e_value_list.append([
                            parameter_dic[query_def][1][h_id]['e-value'],
                            parameter_dic[query_def][1][h_id]['max_score'],
                            parameter_dic[query_def][1][h_id]['identity_of_max_score'],
                            h_id,
                            parameter_dic[query_def][1][h_id]['hit_range_merged']])
                    hit_e_value_list.sort()
        
        
                    hit_garbage_list=[]
                    set_e_value='1e-%s' % str(option_dict["-lim_evalue"])
                    for key in hit_e_value_list:
                        if float(key[0]) > float(set_e_value):
                            hit_garbage_list.append(key[3])
                        if float(key[2]) < float(option_dict['-lim_identity']):
                            hit_garbage_list.append(key[3])
                        if float(key[1]) <float(option_dict['-lim_score']):
                            hit_garbage_list.append(key[3])
                            
                        
                    new_hit_garbage_list=list(set(hit_garbage_list))
                    if len(hit_garbage_list) >= 1:
                        for key in new_hit_garbage_list:
                            del parameter_dic[query_def][1][key]  #hit grabage deletion by e-value
                    new_hit_garbage_list=[] # to clear memory
                    hit_garbage_list=[] # to clear memory
                    hit_e_value_list=[] # to clear memory


#_________________________csv, sequence_dic writing_________________________________________________    


                   #parameter_dic={query_def:'[query_len,{hit_id: {'Hsp_all':hsp_dic,hit_range_merged, hit_def, hit_accession, hit_length,total_score, e-value, mean_identity, query_coverage, subject_coverage}]}'}'}


                    No_result_num=0
    
                    hit_id=list(parameter_dic[query_def][1].keys())
                    evalue_hit_list=[]
                    if len(list(parameter_dic[query_def][1].keys()))==0:
                        No_result_num=No_result_num+1
                        csv_file.writerow([query_def,parameter_dic[query_def][0], "-", "-", "-", "-", "-", "-", "-","-", "-", "-", "-", "-"])

                    temp_dic={}
                    #seq_file=shelve.open(option_dict["-out"]+"_sequence_db.dic") # DB construction for matched query region 
                        
                    if len(list(parameter_dic[query_def][1].keys()))==1:
                        for sub_name in hit_id:
                            evalue_hit_list.append([parameter_dic[query_def][1][sub_name]['e-value'],sub_name])
                        min_hit_list=min(evalue_hit_list)
                        min_hit_id=min_hit_list[1]


                        
                        if len(parameter_dic[query_def][1][min_hit_id]['Hsp_all'])==1:

                            sub_key=list(parameter_dic[query_def][1][min_hit_id]['Hsp_all'].keys())[0]

                            csv_value_list=[
                                query_def,
                                parameter_dic[query_def][0],
                                parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_query-from'],
                                parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_query-to'],
                                parameter_dic[query_def][1][min_hit_id]['hit_def'],
                                parameter_dic[query_def][1][min_hit_id]['hit_length'],
                                parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_hit-from'],
                                parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_hit-to'],
                                parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_bit-score'],
                                parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_evalue'],
                                float(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_identity'])/float(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_align-len']),
                                parameter_dic[query_def][1][min_hit_id]['query_coverage'],
                                parameter_dic[query_def][1][min_hit_id]['subject_coverage']
                                ]
                            if int(option_dict['-seq_parse'])==1:
                                csv_value_list.append(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_qseq'])
                                csv_value_list.append(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_hseq'])

                            csv_file.writerow(csv_value_list)



                        if len(parameter_dic[query_def][1][min_hit_id]['Hsp_all']) > 1:


                            Hsp_id_temp=[]    # To sort hsps by position-start
                            
                            Hsp_keys=list(parameter_dic[query_def][1][min_hit_id]['Hsp_all'].keys())
                            
                            for key in Hsp_keys:
                                Hsp_id_temp.append((parameter_dic[query_def][1][min_hit_id]['Hsp_all'][key]['Hsp_query-from'],key))

                            Hsp_id_temp.sort()

                            Hsp_id=[]    
                            for id in Hsp_id_temp:   
                                Hsp_id.append(id[1])  # To sort hsps by position-end

                            Hsp_id_temp=[]    #To clear memory
                            Hsp_keys=[]       #To clear memory
                            
                            m=0
                            for hsp in Hsp_id:
                                m=m+1
                                csv_value_list=[
                                    query_def,
                                    parameter_dic[query_def][0],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_query-from'],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_query-to'],
                                    parameter_dic[query_def][1][min_hit_id]['hit_def']+'_continued-'+str(m),
                                    parameter_dic[query_def][1][min_hit_id]['hit_length'],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_hit-from'],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_hit-to'],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_bit-score'],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_evalue'],
                                    float(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_identity'])/float(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_align-len']),
                                    parameter_dic[query_def][1][min_hit_id]['query_coverage'],
                                    parameter_dic[query_def][1][min_hit_id]['subject_coverage']
                                    ]

                                if int(option_dict['-seq_parse'])==1:
                                    csv_value_list.append(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_qseq'])
                                    csv_value_list.append(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_hseq'])

                                csv_file.writerow(csv_value_list)

                            
                        

                        


                        #temp_dic[">"+query_def]=parameter_dic[query_def][1][min_hit_id]['Hsp_qseq']# DB writing for matched query region (left)
                        
                        #seq_name_list.write(">"+query_def) # Matched sequence name writing
                        #seq_name_list.write("\n")
                        
                        #temp_dic[">"+parameter_dic[query_def][1][min_hit_id]['hit_def']]=parameter_dic[query_def][1][min_hit_id]['Hsp_hseq'] # DB writing for matched query region (right)
                        
                        #seq_name_list.write(">"+parameter_dic[query_def][1][min_hit_id]['hit_def']) # Matched sequence name writing
                        #seq_name_list.write("\n")

                        #seq_file[">"+query_def]=temp_dic # DB writing for matched query region by dictionary
                        



                    if len(list(parameter_dic[query_def][1].keys())) >= 2:
                        for sub_name in hit_id:
                            evalue_hit_list.append([
                                parameter_dic[query_def][1][sub_name]['e-value'],(1000/parameter_dic[query_def][1][sub_name]['total_score']),parameter_dic[query_def][1][sub_name]['mean_identity'],sub_name])
                        evalue_hit_list.sort()


                        if len(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'])==1:

                            sub_key=list(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'].keys())[0]
                            
                            csv_value_list=[
                                query_def+'_continued-1',
                                parameter_dic[query_def][0],
                                parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_query-from'],
                                parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_query-to'],
                                parameter_dic[query_def][1][evalue_hit_list[0][3]]['hit_def'],
                                parameter_dic[query_def][1][evalue_hit_list[0][3]]['hit_length'],
                                parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_hit-from'],
                                parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_hit-to'],
                                parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_bit-score'],
                                parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_evalue'],
                                float(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_identity'])/float(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_align-len']),
                                parameter_dic[query_def][1][evalue_hit_list[0][3]]['query_coverage'],
                                parameter_dic[query_def][1][evalue_hit_list[0][3]]['subject_coverage']
                                ]
                            if int(option_dict['-seq_parse'])==1:
                                csv_value_list.append(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_qseq'])
                                csv_value_list.append(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_hseq'])

                            csv_file.writerow(csv_value_list)



                        if len(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all']) > 1:

                            Hsp_id_temp=[]    # To sort hsps by position-start
                            Hsp_keys=list(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'].keys())    
                            for key in Hsp_keys:
                                Hsp_id_temp.append((parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][key]['Hsp_query-from'],key))

                            Hsp_id_temp.sort()
                                                                                                                       
                            Hsp_id=[]    
                            for id in Hsp_id_temp:   
                                Hsp_id.append(id[1])  # To sort hsps by position-end

                            Hsp_id_temp=[] # To clear memory
                            Hsp_keys=[]    # To clear memory

                            m=0
                            for hsp in Hsp_id:
                                m=m+1
                                csv_value_list=[
                                    query_def+'_continued-1',
                                    parameter_dic[query_def][0],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_query-from'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_query-to'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['hit_def']+'_continued-'+str(m),
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['hit_length'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_hit-from'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_hit-to'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_bit-score'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_evalue'],
                                    float(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_identity'])/float(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_align-len']),
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['query_coverage'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['subject_coverage']
                                    ]

                                if int(option_dict['-seq_parse'])==1:
                                    csv_value_list.append(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_qseq'])
                                    csv_value_list.append(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_hseq'])

                                csv_file.writerow(csv_value_list)





                        l=0
                        for key in evalue_hit_list[1:]:
                            l=l+1


                            if len(parameter_dic[query_def][1][key[3]]['Hsp_all'])==1:

                                sub_key=list(parameter_dic[query_def][1][key[3]]['Hsp_all'].keys())[0]
                            
                                csv_value_list=[
                                    '%s_continued-%s' % (query_def,str(l+1)),
                                    parameter_dic[query_def][0],
                                    parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_query-from'],
                                    parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_query-to'],
                                    parameter_dic[query_def][1][key[3]]['hit_def'],
                                    parameter_dic[query_def][1][key[3]]['hit_length'],
                                    parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_hit-from'],
                                    parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_hit-to'],
                                    parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_bit-score'],
                                    parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_evalue'],
                                    float(parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_identity'])/float(parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_align-len']),
                                    parameter_dic[query_def][1][key[3]]['query_coverage'],
                                    parameter_dic[query_def][1][key[3]]['subject_coverage']
                                    ]
                                if int(option_dict['-seq_parse'])==1:
                                    csv_value_list.append(parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_qseq'])
                                    csv_value_list.append(parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_hseq'])

                                csv_file.writerow(csv_value_list)



                            if len(parameter_dic[query_def][1][key[3]]['Hsp_all']) > 1:

                                Hsp_id_temp=[]    # To sort hsps by position-start
                                Hsp_keys=list(parameter_dic[query_def][1][key[3]]['Hsp_all'].keys())
                                for key0 in Hsp_keys:
                                    Hsp_id_temp.append((parameter_dic[query_def][1][key[3]]['Hsp_all'][key0]['Hsp_query-from'],key0))
                                Hsp_id_temp.sort()

                                Hsp_id=[]    
                                for id in Hsp_id_temp:   
                                    Hsp_id.append(id[1])  # To sort hsps by position-end

                                Hsp_id_temp=[] #To clear memory
                                Hsp_keys=[]    #To clear memory

                                
                                m=0
                                for hsp in Hsp_id:
                                    m=m+1
                                    csv_value_list=[
                                        '%s_continued-%s' % (query_def,str(l+1)),
                                        parameter_dic[query_def][0],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_query-from'],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_query-to'],
                                        parameter_dic[query_def][1][key[3]]['hit_def']+'_continued-'+str(m),
                                        parameter_dic[query_def][1][key[3]]['hit_length'],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_hit-from'],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_hit-to'],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_bit-score'],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_evalue'],
                                        float(parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_identity'])/float(parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_align-len']),
                                        parameter_dic[query_def][1][key[3]]['query_coverage'],
                                        parameter_dic[query_def][1][key[3]]['subject_coverage']
                                        ]

                                    if int(option_dict['-seq_parse'])==1:
                                        csv_value_list.append(parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_qseq'])
                                        csv_value_list.append(parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_hseq'])

                                    csv_file.writerow(csv_value_list)


                parameter_dic={} # To clear memory
                hsp_dic={} # To clear memory
                hit_parameter_dic={} # To clear memory
                hsp_sub_dic={} # To clear memory
                break









#____________________________Hit Circuit____________________________________
            
            if line.startswith("<Hit>"):

                
                hsp_dic={}
                
                while 1:
                    
                    line=input_file.readline().strip()
                                            
                    if line.startswith("<Hit_id>"):
                        hit_id=line[8:-9]
                        #hit_dic[hit_id]=[] #hit_ID
                    if line.startswith("<Hit_def>"):
                        hit_def=line[9:-10]
                    if line.startswith("<Hit_accession>"):
                        hit_accession=line[15:-16]
                        
                    if line.startswith("<Hit_len>"):
                        hit_length=int(line[9:-10])#hit_length

                             

#________________________________Garbage hsp deletion________________________________________________

                            
                    if line.startswith("</Hit_hsps>"):
                        hsp_dic_keys=list(hsp_dic.keys())
                                         
                        

                        if len(hsp_dic_keys) ==1:
                            set_value=list(range(hsp_dic[hsp_dic_keys[0]]['Hsp_query-from'], hsp_dic[hsp_dic_keys[0]]['Hsp_query-to']+1)) 
                            hsp_hit_merge=set(list(set_value))
                            
                        
                        if len(hsp_dic_keys) > 1:
                            for k in hsp_dic_keys:
                                hsp_dic[k]['hit_range'] = set(list(range(hsp_dic[k]['Hsp_query-from'], hsp_dic[k]['Hsp_query-to']+1))) # 'hit_range' entered to hsp value
                                
                            hsp_score_sort=[]
                            for key in hsp_dic_keys:
                                hsp_score_sort.append([hsp_dic[key]['Hsp_bit-score'], key]) #To sort by score
                            
                            hsp_score_sort.sort()
                            hsp_score_sort.reverse()                                        #To sort by score

                            n=0                                                         #To sort by score
                            hsp_garbage_list=[]                                         #To sort by score
                            for key in hsp_score_sort:                                  #To sort by score
                                n=n+1
                                if n == len(hsp_score_sort):                                #To sort by score
                                    break
                                for key_2 in hsp_score_sort[n:]:
                                    if len(hsp_dic[key[1]]['hit_range'].intersection(hsp_dic[key_2[1]]['hit_range'])) >= int(option_dict['-ov_len']):
                                        hsp_garbage_list.append(key_2[1])

                            temp1=set(hsp_garbage_list)
                            hsp_garbage_new=list(temp1)
                            
                                        
                            if len(hsp_garbage_new) >= 1:
                                for key in hsp_garbage_new:
                                    del hsp_dic[key]  #Hsp garbage deletion by overlap

                            hsp_garbage_list=[] #To clear memory
                            hsp_garbage_new=[] #To clear memory
                            temp1=[] #To clear memory
                            hsp_score_sort=[] #To clear memory

                            ##___________garbage hsp deletion by score_start______________________

                            hsp_dic_keys2=list(hsp_dic.keys())
            
                            if len(hsp_dic_keys2)>1:
                                hsp_score_sort=[]
                                for key in hsp_dic_keys2:
                                    hsp_score_sort.append([hsp_dic[key]['Hsp_bit-score'], float(hsp_dic[key]['Hsp_identity'])/float(hsp_dic[key]['Hsp_align-len']), key])
                                hsp_score_sort.sort()
                                hsp_score_sort.reverse()
                                    
                                hsp_sort_list=hsp_score_sort[1:]
                                hsp_garbage_list=[]
                                hsp_num=1
                                for key in hsp_sort_list:
                                    hsp_num=hsp_num+1
                                    if hsp_num > int(option_dict['-num_hsp']):
                                        hsp_garbage_list.append(key[2])
                                    else:
                                        if float(key[0]) < float(option_dict['-lim_score']):
                                            hsp_garbage_list.append(key[2])
    
                                        if float(key[1]) < float(option_dict['-hsp_identity']):##############################################################
                                            hsp_garbage_list.append(key[2])

                                hsp_garbage_list=list(set(hsp_garbage_list))
                                    
                                if len(hsp_garbage_list) >= 1:
                                    for key in hsp_garbage_list:
                                        del hsp_dic[key]
                                hsp_score_sort=[] #To clear memory
                                hsp_garbage_list=[] #To clear memory
                                hsp_sort_list=[] #To clear memory
                            ##___________garbage hsp deletion by score_end______________________
                                    



                            new_key=list(hsp_dic.keys())

                            if len(new_key) > 1:
                                hsp_hit_merge=set([])
                                for key in new_key:
                                    hsp_hit_merge=set(hsp_dic[key]['hit_range']).union(hsp_hit_merge)
                                
                            if len(new_key) == 1:
                                hsp_hit_merge=hsp_dic[new_key[0]]['hit_range']   #Hsp range merge
                                
                                
                            if len(new_key) == 0:
                                hsp_hit_merge=set([])
                                
                         
                         
                           
                            #hsp_dic={'1':
                                #{'hit_range': "234,2345,34,5346,34,6...",
                                #'Hsp_bit-score': 234235,
                                #'Hsp_query-to':8140,
                                #'Hsp_align-len':3019,
                                #'Hsp_identity':3019,
                                #'Hsp_query-from':5122,
                                #'Hsp_evalue':0.0, 
                                #'Hsp_hit-from':9212,
                                #'Hsp_hseq':'atatgc...',
                                #'Hsp_qseq':'atgc...',
                                #'Hsp_hit-to':12230}}

#____________________________________Value calculate and add___________________________________                        

                        
                        
                        hit_parameter_dic={}
                        hit_parameter_dic['hit_range_merged']=hsp_hit_merge  #merged hit range
                        hit_parameter_dic['hit_def']=hit_def #hit_def value
                        hit_parameter_dic['hit_accession']=hit_accession #hit_accession value
                        hit_parameter_dic['hit_length']=hit_length #hit_length value


                        hsp_keys=list(hsp_dic.keys())
                        
                        for k in hsp_keys:
                            hit_parameter_dic['query_from']=hsp_dic[k]['Hsp_query-from']#########
                            hit_parameter_dic['query_to']=hsp_dic[k]['Hsp_query-to']############
                            hit_parameter_dic['hit_from']=hsp_dic[k]['Hsp_hit-from']############
                            hit_parameter_dic['hit_to']=hsp_dic[k]['Hsp_hit-to']############

                        
                        evalue_list=[]
                        for k in hsp_keys:
                            evalue_list.append(hsp_dic[k]['Hsp_evalue'])
                        hit_parameter_dic['e-value']=min(evalue_list)#e-value
                            
                        
                        hsp_identity=0
                        hsp_length=0
                        hsp_identity_list=[]
                        for k in hsp_keys:
                            hsp_identity=hsp_identity + hsp_dic[k]['Hsp_identity']############
                            hsp_length=hsp_length + hsp_dic[k]['Hsp_align-len']############
                            individual_identity=float(hsp_dic[k]['Hsp_identity'])/float(hsp_dic[k]['Hsp_align-len'])
                            hsp_identity_list.append(individual_identity)
                            
                        hit_parameter_dic['mean_identity']=float(hsp_identity) / float(hsp_length) #mean_identity value############
                        hit_parameter_dic['max_identity']=max(hsp_identity_list) #max_identity value############
                        hit_parameter_dic['min_identity']=min(hsp_identity_list) #min_identity value############

                        total_score=0
                        max_score_list=[]
                        for k in hsp_keys:
                            total_score=total_score + hsp_dic[k]['Hsp_bit-score']
                            max_score_list.append([float(hsp_dic[k]['Hsp_bit-score']),float(hsp_dic[k]['Hsp_identity'])/float(hsp_dic[k]['Hsp_align-len'])])
                        hit_parameter_dic['total_score']=total_score #total_score value
                        hit_parameter_dic['max_score']=max(max_score_list)[0]#max_score value
                        hit_parameter_dic['identity_of_max_score']=max(max_score_list)[1]#identity of max_score value


                        red_len_query=0
                        split_length=[]
                        for k in hsp_keys:
                            hit_covered_len_query = hsp_dic[k]['Hsp_query-to'] - hsp_dic[k]['Hsp_query-from'] + 1
                            
                            hit_parameter_dic['query_coverage']=int(hit_covered_len_query)  #query_coverage value

                        hit_covered_len_subject=0
                        for k in hsp_keys:
                            hit_covered_len_subject = hit_covered_len_subject + (hsp_dic[k]['Hsp_hit-to'] - hsp_dic[k]['Hsp_hit-from'] + 1)
                        hit_parameter_dic['subject_coverage'] = abs(int(hit_covered_len_subject)) #subject_coverage value

                        hit_parameter_dic['Hsp_qseq']=hsp_dic[k]['Hsp_qseq'] # Aligned query sequence
                        hit_parameter_dic['Hsp_hseq']=hsp_dic[k]['Hsp_hseq'] # Aligned hit sequence
                        hit_parameter_dic['Hsp_all']=hsp_dic

    
                        parameter_dic[query_def][1][hit_id]=hit_parameter_dic #{query_def:'[query_len,{hit_id: {'Hsp_all':hsp_dic,hit_range_merged, hit_def, hit_accession, hit_length,total_score, e-value, mean_identity, query_coverage, subject_coverage}]}'}'}


                        break
                

#_____________________________________Hsp Circuit__________________________

                    if line.startswith("<Hsp>"):
                        hsp_sub_dic={} #{hit_range, hsp_bit-score,hsp_evalue,hsp_query-from,hsp_query-to,hsp_hit-from,hsp_hit-to,hsp_identity, hsp_align-len}
                        while 1:
                            line=input_file.readline().strip()
                            if line.startswith("<Hsp_num>"):
                                hsp_num=line[9:-10]
                            if line.startswith("<Hsp_bit-score>"):
                                hsp_sub_dic['Hsp_bit-score']=float(line[15:-16])
                            if line.startswith("<Hsp_evalue>"):
                                hsp_sub_dic['Hsp_evalue']=float(line[12:-13])
                            if line.startswith("<Hsp_query-from>"):
                                hsp_sub_dic['Hsp_query-from']=int(line[16:-17])
                            if line.startswith("<Hsp_query-to>"):
                                hsp_sub_dic['Hsp_query-to']=int(line[14:-15])
                            if line.startswith("<Hsp_hit-from>"):
                                hsp_sub_dic['Hsp_hit-from']=int(line[14:-15])
                            if line.startswith("<Hsp_hit-to>"):
                                hsp_sub_dic['Hsp_hit-to']=int(line[12:-13])
                            if line.startswith("<Hsp_identity>"):
                                hsp_sub_dic['Hsp_identity']=int(line[14:-15])
                            if line.startswith("<Hsp_align-len>"):
                                hsp_sub_dic['Hsp_align-len']=int(line[15:-16])
                            if line.startswith("<Hsp_qseq>"):
                                Hsp_qseq_list=[]
                                while 1:
                                    Hsp_qseq_list.append(line)
                                    line=input_file.readline().strip()
                                    if line.startswith("<Hsp_hseq>"):
                                        Hsp_qseq="".join(Hsp_qseq_list)[10:-11]#.replace("-","") # Conserve the alignment format
                                        hsp_sub_dic['Hsp_qseq']=Hsp_qseq
                                        break
                            if line.startswith("<Hsp_hseq>"):
                                Hsp_hseq_list=[]
                                while 1:
                                    Hsp_hseq_list.append(line)
                                    line=input_file.readline().strip()
                                    if line.startswith("<Hsp_midline>"):
                                        Hsp_hseq="".join(Hsp_hseq_list)[10:-11]#.replace("-","") # Conserve the alignment format
                                        hsp_sub_dic['Hsp_hseq']=Hsp_hseq
                                        break


                            if line.startswith("</Hsp>"):
                                hsp_dic[hsp_num]=hsp_sub_dic
                                break
                        







                                 

