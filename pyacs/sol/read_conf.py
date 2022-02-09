
class Conf:
    
    """Configuration: various parameters read by the program"""

    
    def __init__ (self, conf_file=None,verbose=False):

        import logging
        import pyacs.message.message as MESSAGE
        import pyacs.message.verbose_message as VERBOSE
        import pyacs.message.error as ERROR
        import pyacs.message.warning as WARNING
        import pyacs.message.debug_message as DEBUG

        VERBOSE("Initializing default values for PYACS")
        self.ref_exclude_all=[]
        self.ref_reject={}
        self.rename={}
        self.ref_only_all = []
        self.min_repeat=[]
        
        self.ref_outlier_up_value=50.0 # 50 mm for up outlier detection 
        self.ref_outlier_en_value=30.0 # 30 mm for e&n outlier detection 
        self.ref_prefit_value=10.0     # 10 m 
        if conf_file:
            MESSAGE("Reading configuration file: %s" % conf_file)
            fs = open(conf_file, 'r')
            for line in fs:
                if (len(line)<2 or line[0] == '#' or line[0] == '*'):continue
                else:
                    lline=line.split()

                    # min_repeats
                    if (lline[0] == 'min_repeat'):
                        lline.pop(0)
                        self.min_repeat.extend(lline) 

                    # ref_exclude_all
                    if (lline[0] == 'ref_exclude_all'):
                        lline.pop(0)
                        self.ref_exclude_all.extend(lline) 
                    
                    # ref_only_all
                    if (lline[0] == 'ref_only_all'):
                        lline.pop(0)
                        self.ref_only_all.extend(lline) 
                    
                    # ref_outlier_up_value
                    if (lline[0] == 'ref_outlier_up_value'):
                        self.ref_outlier_up_value=float(lline.pop()) 
                    
                    # ref_outlier_en_value
                    if (lline[0] == 'ref_outlier_en_value'):
                        self.ref_outlier_en_value=float(lline.pop()) 
                    
                    # ref_prefit_value
                    if (lline[0] == 'ref_prefit_value'):
                        self.ref_prefit_value=float(lline.pop()) 
                    
                    # rename
                    if (lline[0] == 'rename'):
                        code_orig=lline[1]
                        code_new=lline[2]
                        if (len(lline)>3):
                            for i in range(3,len(lline)):
                                sinex=lline[i]
                                if sinex not in list(self.rename.keys()):self.rename[sinex]=[]
                                self.rename[sinex].append((code_orig,code_new))
                        else:
                            if not 'all' in list(self.rename.keys()):self.rename['all']=[]
                            self.rename['all'].append((code_orig,code_new))

                    # reject_ref
                    if (lline[0] == 'reject_ref'):
                        code=lline[1]
                        component=lline[2]
                        sinex=lline[3]
                        if sinex not in list(self.ref_reject.keys()): 
                            self.ref_reject[sinex]={}
                            self.ref_reject[sinex][code]=component
                        else:
                            if code not in self.ref_reject[sinex]:
                                self.ref_reject[sinex][code]=component
                            else:
                                if component not in self.ref_reject[sinex][code]:
                                    self.ref_reject[sinex][code]=self.ref_reject[sinex][code]+component
        else:
            MESSAGE("No configuration file provided; Using defaults values")

def main():
    
    conf=Conf()
    for k, v in conf.__dict__.items():
        if v==[] or v=={}:
            print("#",k, v)
        else:
            print(k,v)

if __name__ == '__main__':
    main()

