import random
prims1 = gfx.PrimList("prims one")
prims2 = gfx.PrimList("prims two")
for i in range(1000):
    prims1.AddLine(10.0*geom.Vec3(random.random()-1,
                                 random.random()-1,
                                 random.random()-0.5),
                  10.0*geom.Vec3(random.random()-0.5,
                                 random.random()-0.5,
                                 random.random()-0.5),
                  gfx.Color(random.random(),
                            random.random(),
                            0.5)
                  );
    prims2.AddLine(10.0*geom.Vec3(random.random(),
                                 random.random(),
                                 random.random()-0.5),
                  10.0*geom.Vec3(random.random()-0.5,
                                 random.random()-0.5,
                                 random.random()-0.5),
                  gfx.Color(random.random(),
                            0.5,
                            random.random())
                  );
scene.Add(prims1)
scene.Add(prims2)
scene.SetCenter(geom.Vec3(0,0,0))



                            
                  
                            
