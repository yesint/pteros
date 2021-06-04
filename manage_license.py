import os,sys,re

rootdir = '.'

license = open("license_header.txt").read()

for root, subdirs, files in os.walk(rootdir):
    for f in files:
        if f.endswith(".cpp") or f.endswith(".h"):
            # Try this file
            try:
                text = open(os.path.join(root,f)).read()
            except:
                pass
            # Try to extract leading comment
            lines = text.split('\n');

            if not lines[0].startswith("/*"):
                continue

            body = ""
            good = False
            in_body = False
            for line in lines:
                if not good and not in_body and re.search("Pteros",line):
                    good = True
                
                if not in_body and line.endswith("*/"):
                    in_body = True
                    continue
                    
                if good and in_body:
                    body += line+"\n"
    
            if good:
                print(os.path.join(root,f))
                
                with open(os.path.join(root,f),'w') as out:
                    out.write(license+body)
        
