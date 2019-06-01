from pyrosetta_helper import pose_structure_df

class BluePrintEntity():
    '''
     1   V  LE  R   (position 1, Val, loop, abego type E, rebuild)
     0   V  EX  R   (insert a residue, Val, sheet, any abego, rebuild)
     2   V  EB  .   (position 2, Val, sheet, abego type B, do not rebuild)
    '''
    def __init__(self,position,aa,ss,abego,rebuild):
        self.position = position
        self.aa = aa
        self.ss = ss
        self.abego = abego
        self.rebuild = rebuild
       
    def set_position(self,p):
        self.position = p
    
    def set_aa(self,aa):
        self.aa = aa
        
    def set_ss(self,ss):
        self.ss = ss
        
    def set_abego(self,abego):
        self.abego = abego
        
    def set_rebuild(self,r):
        self.rebuild = r
        
    def get_line(self):
        return "{position} {aa} {ss}{abego} {rebuild}".format(
            position=self.position,
            aa=self.aa,
            ss=self.ss,
            abego=self.abego,
            rebuild='R' if self.rebuild else '.')
    def __str__(self):
        return "<"+self.get_line()+">"
    
    def __repr__(self):
        return self.__str__()
        
class BluePrintManager():
    
    def __init__(self,pose=''):
    
        self.blueprint = {}
        self.blueprint_size = 0
        if pose:
            self.pose = pose
            self.blueprint_size = self.pose.size()
            self._make_blueprint_from_pose()
            
    def _make_blueprint_from_pose(self):        
        pose_df = pose_structure_df(self.pose)
        index = 1
        for key,val in pose_df.iterrows():
            #print(key,val)
            position = key
            aa = val['Residue']
            abego = val['ABEGO']
            ss = val['SS']
            rebuild = False
            entity = BluePrintEntity(position,aa,ss,abego,rebuild)
            self.blueprint[index] = entity
            index += 1
            
    def set_rebuild(self,index,rebuild):
        self.blueprint[index].set_rebuild(rebuild)
        
    def set_aa(self,index,aa):
        self.blueprint[index].set_aa(aa)
    
    def set_ss(self,index,ss):
        self.blueprint[index].set_ss(ss)
        
    def set_abego(self,index,abego):
        self.blueprint[index].set_abego(abego)
    
    def insert_aa_after(self,after,aa,ss,abego):
        old_blueprint = self.blueprint.copy()
        self.blueprint = {}
        entity = BluePrintEntity(0,aa,ss,abego,True)
        for index in old_blueprint:
            #pdb.set_trace()
            if index < after:
                self.blueprint[index] = old_blueprint[index]
            elif index == after:
                self.blueprint[index]=old_blueprint[index]

                self.blueprint[index+1] = entity
            else:
                self.blueprint[index+1]=old_blueprint[index]
                
    def insert_aa_stretch_after(self,after,stretch):
        
        for aa in stretch[::-1]:
            self.insert_aa_after(after,aa,'L','X')
            
    def insert_aa_before(self,before,aa,ss,abego):
        old_blueprint = self.blueprint.copy()
        self.blueprint = {}
        entity = BluePrintEntity(0,aa,ss,abego,True)
        for index in old_blueprint: 
            #pdb.set_trace()
            if index > before:
                self.blueprint[index+1] = old_blueprint[index]
                continue
            elif index == before: 
                self.blueprint[index] = entity
                self.blueprint[index+1]=old_blueprint[index]
            
            else:
                self.blueprint[index] = old_blueprint[index]

                
    def insert_aa_stretch_before(self,before,stretch):
        for aa in stretch[::-1]:
            self.insert_aa_before(before,aa,'L','X')
            
    def writefile(self,filename):
        with open(filename,'w') as f:
            f.write(self.__str__())
    
    def __str__(self):
       return '\n'.join([i.get_line() for i in self.blueprint.values()])
        
        
    def __repr__(self):
        return "BLUEPRINT:\n"+'\n'.join(
            ["Index={} Line={}".format(key,b.__str__()) for key,b in self.blueprint.items()])