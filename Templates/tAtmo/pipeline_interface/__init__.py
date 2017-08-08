__DESCRIPTION__=""" exoplanet pipeline par file parser
By M.Maris - 0.0 - 2017 June 27 -
"""
class Singleton(object):
  """
  Class implementing the Singleton design pattern. There are several ways to implement
  the singleton pattern in python. The most elegant way is to use the decorators but such
  construct has been introduced only starting from python 2.6
  """
  def __new__(cls, *args, **kwds):
    if not '_the_instance' in cls.__dict__:
      cls._the_instance =  object.__new__(cls)
      cls._the_instance.init(*args, **kwds)
    return cls._the_instance
  def init(self, *args, **kwds):
    pass

class _config_setup(Singleton) :
   "RELEASE related information"
   def init(self) :
      pass
   def __call__(self,inputString,isHereDoc=False) :
      import time
      self._instantiated_at_=time.asctime()
      self._keywords_=[]
      for k in (\
               open(inputString,'r') \
            if isHereDoc==False else (\
                  inputString.split('\n') 
               if type(inputString) == type('') else \
                  inputString \
                  )\
               ):
         l=k.strip()
         if len(l) > 0 :
            if l[0]!='#' :
               ll=l.split('=')
               if ll[1].strip()!='' :
                  self[ll[0].strip()]=ll[1].strip() 
                  self._keywords_.append(ll[0].strip())
   #
   def copy(self) :
      import copy
      return copy.deepcopy(self)
   #
   def keys(self) :
      return self._keywords_
   #
   def __setitem__(self,this,that) :
      self.__dict__[this]=that
   #
   def __getitem__(self,this) : 
      if len(self.__dict__) == 0 : return
      return self.__dict__[this]
   #
   def __str__(self) :
      l=[]
      for k in self.keys() :
         l.append('%s = %s'%(k,self[k]))
      return "\n".join(l)
   #
   def top_pressure(self) :
      if not "ptop" in self.keys() : return 0.01
      return float(self['ptop'].strip())
   #
   def OverWriteParFile(self) :
      if not "overWriteParFile" in self.keys() : return False
      ss = str(self["overWriteParFile"]).strip().lower()
      if len(ss) == 0 : return False
      if ss[0] == '1' or ss[0]=='t' or ss[0]=='true' : return True
      return False
   #
   def banner(self) :
      print 
      print "***"
      print '***instantiated_at : ',self._instantiated_at_
      print "***"
      print self
      print
CONFIG=_config_setup()

class exop_par_file :
   def __init__(self,*arg) :
      import numpy as np
      self.fname=None
      self.nlines=0
      self.nline=[]
      self.line=[]
      self.key=[]
      self.value=[]
      self.comment=[]
      self.type=[]
      self.iscomment=[]
      if len(arg) == 1 : 
         self.fname=arg[0]
         for _k in open(self.fname,'r') :
            k=_k.split('\n')[0]
            self.line.append(k)
            self.nline.append(self.nlines)
            self.nlines+=1
            if k.strip()=='' :
               kk=['','','','']
            else :
               kk=k.split('!')
               if len(kk) == 2 :
                  kk.append('')
                  kk.append('')
            self.iscomment.append(kk[0].strip().lower()=='comment')
            self.key.append(kk[0].strip())
            self.value.append(kk[1].strip())
            self.comment.append(kk[2].strip())
            self.type.append('' if kk[0].strip().lower()=='comment' else kk[3].strip())
      for k in ['nline','line','key','value','comment','type','iscomment'] :
         self.__dict__[k]=np.array(self.__dict__[k])
   #      
   def __len__(self) : return self.nlines
   #      
   def copy(self) :
      import copy
      return copy.deepcopy(self)
   #      
   def keys(self) : return self.key
   #      
   def argselect(self,this) :
      import numpy as np
      idx= np.where(self.key==this)[0]
      if len(idx) == 0 : return -1
      return idx[0]
   #      
   def __getitem__(self,this) :
      import numpy as np
      idx=self.argselect(this)
      return None if idx < 0 else (self.nline[idx],self.key[idx],self.value[idx],self.comment[idx],self.type[idx])
   #
   def update(self,name,value) :
      idx=self.argselect(name)
      if idx < 0 : return 
      self.value[idx]=value
   #
   def __line_compose(self,key,value,comment,Type) :
      import numpy as np
      if key.strip().lower()=='comment' :
         return '!'.join([key,'',comment,''])
      if key.strip()=='' :
         return ''
      return '!'.join([key,value,comment,Type])
   #
   def append(self,key,value,comment,Type) :
      import numpy as np
      if key.strip().lower()=='comment' :
         self.nline = np.concatenate([self.nline,[self.nline[-1]+1]])
         self.key = np.concatenate([self.key,[key.strip().upper()]])
         self.value = np.concatenate([self.value,['']])
         self.comment = np.concatenate([self.comment,[comment]])
         self.type = np.concatenate([self.type,['']])
         self.iscomment = np.concatenate([self.iscomment,[True]])
         self.line = np.concatenate([self.line,[self.__line_compose(key.strip().upper(),'',comment,'')]])
         self.nlines+=1
         return
      idx=np.where(self.iscomment)[0]
      if len(idx) == 0:
         self.nline = np.concatenate([self.nline,[self.nline[-1]+1]])
         self.key = np.concatenate([self.key,[key.strip().upper()]])
         self.value = np.concatenate([self.value,[value]])
         self.comment = np.concatenate([self.comment,[comment]])
         self.type = np.concatenate([self.type,[Type.strip().upper()]])
         self.iscomment = np.concatenate([self.iscomment,[False]])
         self.line = np.concatenate([self.line,[self.__line_compose(key.strip().upper(),value,comment,Type.strip().upper())]])
         self.nlines+=1
      else :
         il=idx[0]-1
         ifc=idx[0]
         #return self.nline[:ifc],[self.nline[ifc]],self.nline[ifc:]
         self.nline = np.concatenate([self.nline[:ifc],[self.nline[ifc]],self.nline[ifc:]+1])
         self.key = np.concatenate([self.key[:ifc],[key.strip().upper()],self.key[ifc:]])
         self.value = np.concatenate([self.value[:ifc],[value],self.value[ifc:]])
         self.comment = np.concatenate([self.comment[:ifc],[comment],self.comment[ifc:]])
         self.type = np.concatenate([self.type[:ifc],[Type.strip().upper()],self.type[ifc:]])
         self.iscomment = np.concatenate([self.iscomment[:ifc],[False],self.iscomment[ifc:]])
         self.line= np.concatenate([self.line[:ifc],[self.__line_compose(key.strip().upper(),value,comment,Type.strip().upper())],self.line[ifc:]])
         self.nlines+=1
   #
   def __str__(self) :
      out=[]
      for k in range(len(self)) :
         line="%4d:%s: "%(self.nline[k],'*' if self.iscomment[k] else ' ')
         line+=self.key[k]
         line+='!'+self.value[k]
         line+='!'+self.comment[k]
         line+='' if self.iscomment[k] else (' !'+self.type[k])
         out.append(line)
      return '\n'.join(out)
   #
   def banner(self) :
      print 
      print EPF.fname
      print 
      print "Number of lines : ",len(self)
      print
      print "Content:"
      print
      print self
      print
   #
   def tofile(self,ofile) :
      oo=open(ofile,'w')
      for k in self.line :
         oo.write(k+'\n')
      oo.close()

if __name__=="__main__" : 
   print 
   print "example of usage"
   print
   #   
   EPF=exop_par_file('exoplanets-test.par')
   EPF.banner()
   #   
   print "Add parameter"
   EPF.append('mio_par','mio_value','my parameter','str')
   #   
   print "Add comment"
   EPF.append('comment','','added my parameter','')
   #   
   EPF.banner()
   #   
   print "Write to file test.par"
   EPF.tofile('test.par')
   #   
   print "Example ended"
   print

   