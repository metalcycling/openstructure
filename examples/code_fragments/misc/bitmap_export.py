# this works with ost_cl, too !

import random
prims = gfx.PrimList("test")
for i in range(1000):
    prims.AddLine(50.0*geom.Vec3(random.random(),
                                 random.random(),
                                 random.random()),
                  50.0*geom.Vec3(random.random(),
                                 random.random(),
                                 random.random()),
                  gfx.Color(random.random(),
                            random.random(),
                            random.random())
                  );

scene.Add(prims)
scene.SetCenter(prims.GetCenter())
scene.Autoslab()
scene.Export("test_bitmap_export.png", 500, 500)


                            
                  
                            
