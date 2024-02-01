import unittest,os
from ost.db import LinearIndexer
from ost.db import LinearCharacterContainer
from ost.db import LinearPositionContainer
from ost import geom


def ComparePosList(l1, l2):
  if len(l1)!= len(l2):
    return False
  eps = 0.0001
  for i in range(len(l1)):
    if geom.Distance(l1[i],l2[i]) > eps:
      return False
  return True


class TestLinearDB(unittest.TestCase):

  def testLinearIndexer(self):

    indexer = LinearIndexer()
    self.assertEqual(indexer.GetNumResidues(), 0)

    # check whether error is thrown when list parameters do not have
    # same length
    self.assertRaises(Exception, indexer.AddAssembly, "HELLOU",
                      ["AA", "B"], [20, 33, 1])

    # correct it and add first assembly
    indexer.AddAssembly("HELLOU", ["AA", "B", "CCC"], [20, 33, 1])

    # check whether according entry has been added
    self.assertEqual(indexer.GetAssemblies(), ["HELLOU"])
    self.assertEqual(indexer.GetNumResidues(), 54)

    # check whether error get thrown when the same assembly is added again
    # and directly add another assembly
    self.assertRaises(Exception, indexer.AddAssembly, "HELLOU",  
                      ["AA", "B", "CCC"], [20, 33, 1])

    indexer.AddAssembly("HELLOUYOU", ["X", "Y", "ZZZ"], [10, 26, 20])

    # check whether according entry has been added. we need to sort the
    # output since the GetAssemblies function gives no guarantee about ordering
    entries = indexer.GetAssemblies()
    entries.sort()
    self.assertEqual(entries, ["HELLOU", "HELLOUYOU"])
    self.assertEqual(indexer.GetNumResidues(), 110)

    # check the GetDataRange functions and all their errors when 
    # inexistent data is accessed
    self.assertRaises(Exception, indexer.GetDataRange, "HELL", "CCC")
    self.assertRaises(Exception, indexer.GetDataRange, "HELLOU", "CC")
    self.assertEqual(indexer.GetDataRange("HELLOU", "CCC"), (53,54))

    self.assertRaises(Exception, indexer.GetDataRange, "HELL")
    self.assertEqual(indexer.GetDataRange("HELLOU"), (0,54))

    # add some more assemblies to test proper removal
    indexer.AddAssembly("THIRDENTRY", ["X"], [5])
    indexer.AddAssembly("FOURTHENTRY", ["R","asdf"],[22, 10])

    # check the errors to be thrown when the RemoveAssembly function is called
    self.assertRaises(Exception, indexer.RemoveAssembly, "HELLOUYO")

    # before removing anything we check tons of ranges, we expect
    # them to be ordered as they have been added to the indexer
    self.assertEqual(indexer.GetDataRange("HELLOU"), (0, 54))
    self.assertEqual(indexer.GetDataRange("HELLOU", "B"), (20, 53))    
    self.assertEqual(indexer.GetDataRange("HELLOUYOU"), (54, 110))
    self.assertEqual(indexer.GetDataRange("HELLOUYOU", "Y"), (64, 90))    
    self.assertEqual(indexer.GetDataRange("THIRDENTRY"), (110, 115))
    self.assertEqual(indexer.GetDataRange("THIRDENTRY", "X"), (110, 115))  
    self.assertEqual(indexer.GetDataRange("FOURTHENTRY"), (115, 147))
    self.assertEqual(indexer.GetDataRange("FOURTHENTRY", "asdf"), (137, 147))  

    self.assertEqual(indexer.GetNumResidues(), 147)    

    # remove an assembly
    indexer.RemoveAssembly("HELLOUYOU")

    # we expect, that the data ranges for assembly HELLOU stays the same,
    # but for all others there should be a decrease of 10+26+20=56
    self.assertEqual(indexer.GetDataRange("HELLOU"), (0, 54))
    self.assertEqual(indexer.GetDataRange("HELLOU", "B"), (20, 53))      
    self.assertEqual(indexer.GetDataRange("THIRDENTRY"), (110 - 56, 115 - 56))
    self.assertEqual(indexer.GetDataRange("THIRDENTRY", "X"), (110 - 56, 115 - 56))  
    self.assertEqual(indexer.GetDataRange("FOURTHENTRY"), (115 - 56, 147 - 56))
    self.assertEqual(indexer.GetDataRange("FOURTHENTRY", "asdf"), (137 - 56, 147 - 56))  
    self.assertEqual(indexer.GetNumResidues(), 147 - 56)

    # We additionally check, whether the thing is really gone
    current_entries = indexer.GetAssemblies()
    current_entries.sort()
    self.assertEqual(current_entries, ["FOURTHENTRY", "HELLOU", "THIRDENTRY"])

    # Let's finally test saving and loading
    indexer.Save("test_indexer.dat")
    loaded_indexer = LinearIndexer.Load("test_indexer.dat")
    self.assertEqual(indexer, loaded_indexer)

    if os.path.exists("test_indexer.dat"):
        os.remove("test_indexer.dat")


  def testCharacterContainer(self):

    container = LinearCharacterContainer()

    self.assertEqual(container.GetNumElements(), 0)

    # define a few sequences to play with
    s1 = "AAAAAAA"
    s2 = "BBB"
    s3 = "CCCCCCCCCC"

    # add some stuff
    container.AddCharacters(s1)
    container.AddCharacters(s2)
    container.AddCharacters(s3)

    # check whether stuff has properly been added
    self.assertEqual(container.GetNumElements(), len(s1) + len(s2) + len(s3))
    self.assertEqual(container.GetCharacters((0, container.GetNumElements())), 
                     s1 + s2 + s3)
    self.assertEqual(container.GetCharacters((0, len(s1))), "AAAAAAA")
    self.assertEqual(container.GetCharacters((len(s1), len(s1) + len(s2))), "BBB")
    self.assertEqual(container.GetCharacters((len(s1) + len(s2), 
                                              len(s1) + len(s2) + len(s3))), 
                                              "CCCCCCCCCC")


    # check whether errors get thrown for invalid ranges
    self.assertRaises(Exception, container.GetCharacters, (5,5))
    self.assertRaises(Exception, container.GetCharacters, (5,4))
    self.assertRaises(Exception, container.GetCharacters, (5,100))
    self.assertRaises(Exception, container.GetCharacters, (100,120))
    self.assertRaises(Exception, container.ClearRange, (5,5))
    self.assertRaises(Exception, container.ClearRange, (5,4))
    self.assertRaises(Exception, container.ClearRange, (5,100))
    self.assertRaises(Exception, container.ClearRange, (100,120))

    # let's remove s2 and then check whether the content is as expected
    container.ClearRange((len(s1),len(s1) + len(s2)))
    self.assertEqual(container.GetNumElements(), len(s1) + len(s3))
    self.assertEqual(container.GetCharacters((0, len(s1))), "AAAAAAA")
    self.assertEqual(container.GetCharacters((len(s1), len(s1) + len(s3))), "CCCCCCCCCC")

    # let's check the single character extraction function
    full_expected_string = s1 + s3
    for i in range(container.GetNumElements()):
      self.assertEqual(full_expected_string[i], container.GetCharacter(i))

    # save / load and check for equality
    container.Save("test_character_container.dat")
    loaded_container = LinearCharacterContainer.Load("test_character_container.dat")
    self.assertEqual(container, loaded_container)

    if os.path.exists("test_character_container.dat"):
        os.remove("test_character_container.dat")


  def testPositionContainer(self):

    container = LinearPositionContainer()

    self.assertEqual(container.GetNumElements(), 0)

    # define a few position lists to play with
    pos_1 = geom.Vec3List()
    pos_2 = geom.Vec3List()
    pos_3 = geom.Vec3List()

    pos_1.append(geom.Vec3(-10.234, -45.333, -45.338))
    pos_1.append(geom.Vec3(10.234, 45.333, 45.338))
    pos_1.append(geom.Vec3(10.234, 45.333, 45.338))
    pos_1.append(geom.Vec3(10.234, 45.333, 45.338))
    pos_1.append(geom.Vec3(10.234, 45.333, 45.338))
    pos_1.append(geom.Vec3(10.234, 45.333, 45.338))

    pos_2.append(geom.Vec3(50000.34, 32.43, 33.2222))
    pos_2.append(geom.Vec3(-50000.34, -32.43, -33.2222))
    pos_2.append(geom.Vec3(50000.34, 32.43, 33.2222))

    pos_3.append(geom.Vec3(87.909, 65.222, 1000.556))
    pos_3.append(geom.Vec3(87.909, 65.222, 1000.556))

    # get the rounded target values (the accuracy with which they)
    # are stored in the container)
    pos_1_target = geom.Vec3(10.23, 45.33, 45.34)
    pos_2_target = geom.Vec3(50000.3, 32.4, 33.2) # reduced accuracy!
    pos_3_target = geom.Vec3(87.91, 65.22, 1000.56)

    pos_1_target_list = geom.Vec3List()
    pos_2_target_list = geom.Vec3List()
    pos_3_target_list = geom.Vec3List()

    for i in range(6):
      pos_1_target_list.append(pos_1_target)
    pos_1_target_list[0] = -pos_1_target_list[0]

    for i in range(3):
      pos_2_target_list.append(pos_2_target)
    pos_2_target_list[1] = -pos_2_target_list[1]

    pos_3_target_list.append(pos_3_target)
    pos_3_target_list.append(pos_3_target)

    # add it
    container.AddPositions(pos_1)
    container.AddPositions(pos_2)
    container.AddPositions(pos_3)

    extracted_pos_list = geom.Vec3List()
    extracted_pos = geom.Vec3()

    # check whether stuff has properly been added
    self.assertEqual(container.GetNumElements(), 
                     len(pos_1) + len(pos_2) + len(pos_3))

    container.GetPositions((0,len(pos_1)), extracted_pos_list)

    self.assertTrue(ComparePosList(pos_1_target_list, extracted_pos_list))

    container.GetPositions((len(pos_1),len(pos_1)+len(pos_2)), extracted_pos_list)
    self.assertTrue(ComparePosList(pos_2_target_list, extracted_pos_list))

    container.GetPositions((len(pos_1)+len(pos_2),
                            len(pos_1)+len(pos_2)+len(pos_3)), extracted_pos_list)
    self.assertTrue(ComparePosList(pos_3_target_list, extracted_pos_list))

    # check whether errors get thrown for invalid ranges
    dummy_list = geom.Vec3List()
    dummy_pos = geom.Vec3()
    self.assertRaises(Exception, container.GetPositions, (5,5), dummy_list)
    self.assertRaises(Exception, container.GetPositions, (5,4), dummy_list)
    self.assertRaises(Exception, container.GetPositions, (5,100), dummy_list)
    self.assertRaises(Exception, container.GetPositions, (100,120), dummy_list)
    self.assertRaises(Exception, container.GetPosition, 100, dummy_pos)
    self.assertRaises(Exception, container.ClearRange, (5,5))
    self.assertRaises(Exception, container.ClearRange, (5,4))
    self.assertRaises(Exception, container.ClearRange, (5,100))
    self.assertRaises(Exception, container.ClearRange, (100,120))

    # lets remove pos_2_list and check, whether the result is as expected
    container.ClearRange((len(pos_1), len(pos_1)+len(pos_2)))
    
    self.assertEqual(container.GetNumElements(), len(pos_1) + len(pos_3))

    container.GetPositions((0,len(pos_1)), extracted_pos_list)
    self.assertTrue(ComparePosList(pos_1_target_list, extracted_pos_list))

    container.GetPositions((len(pos_1), len(pos_1) + len(pos_3)), extracted_pos_list)
    self.assertTrue(ComparePosList(pos_3_target_list, extracted_pos_list))

    # check the single position getter
    pos_13_target_list = geom.Vec3List()
    for item in pos_1_target_list:
      pos_13_target_list.append(item)
    for item in pos_3_target_list:
      pos_13_target_list.append(item)

    eps = 0.00001
    for i in range(len(pos_13_target_list)):
      container.GetPosition(i, extracted_pos)
      self.assertTrue(geom.Distance(extracted_pos, pos_13_target_list[i]) < eps)

    # save / load and check for equality
    container.Save("test_position_container.dat")
    loaded_container = LinearPositionContainer.Load("test_position_container.dat")
    self.assertEqual(container, loaded_container)

    if os.path.exists("test_position_container.dat"):
        os.remove("test_position_container.dat")


if __name__ == "__main__":

  from ost import testutils
  testutils.RunTests()

