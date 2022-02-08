# -*- coding: utf-8 -*-

import time
from linkbudgetcalc import *
class EventDetector():
  def __init__(self, eventType, inputList, time=time.ctime(time.time())):
    self.eventType = eventType
    self.time = time
    self.inputList = inputList

class EventLogger():
  # eventType = integer 1-6
  def __init__(self, fileName="EventLogger.txt"):
    self.fileName = fileName
    t = time.ctime(time.time())
    logger = open(fileName,"w")
    logger.write(t + "\tMISSION START\n")
  
  # Log a nominal information about what immediate event happened
  def info(self, obj, msg):
    time = obj.time
    objtype = obj.eventType
    logger = open(self.fileName,"a")
    logger.writelines([time, "\tINFO".ljust(10), objtype, ": ", msg, "\n"])
    logger.close()
    return
  
  # Log the start of a continuous event
  def start(self, obj):
    time = obj.time
    objtype = obj.eventType
    logger = open(self.fileName,"a")
    logger.writelines([time, "\tSTART".ljust(10), objtype, ": ", objtype, " started to launch", "\n"])
    logger.close()
    return

  # Log the end of a continuous event
  def stop(self, obj):
    time = obj.time
    objtype = obj.eventType
    logger = open(self.fileName,"a")
    logger.writelines([time, "\tSTOP".ljust(10), objtype, ": ", objtype, " stopped to launch", "\n"])
    logger.close()
    return

  # Log a off-nominal event happened
  def warning(self, obj, msg):
    time = obj.time
    objtype = obj.eventType
    logger = open(self.fileName,"a")
    logger.writelines([time, "\tWARNING".ljust(10), objtype, ": ", msg, "\n"])
    logger.close()
    return
  
  # Find a event history containing "item", return "not found" if not existed
  def trackHistory(self, item):
    logger = open(self.fileName, "r")

    flag = 0
    num = 0

    for line in logger:
      if item in line:
        flag = 1
        num += 1
        print(line)

    if flag == 0:
      print('Event', '\"' + item + '\"', 'NOT FOUND')
    else:
      print(str(num) + " records found")

    logger.close()
    return

det1 = EventDetector("test", [1, 2, 3])
print(det1.time)
print()
log1 = EventLogger()
det2 = EventDetector("Comms", [])
det3 = EventDetector("Orbit", [])
log1.info(det2,"This comms object was created")
log1.start(det3)
log1.stop(det3)
log1.warning(det2, "Something odd happened")

log1.trackHistory("Orbit")
log1.trackHistory("Test")


import itertools
from itertools import combinations
class EventHandler():
  def __init__(self, eventType, viewLimit, interest, powerLimit):
    self.eventType = eventType
    self.viewLimit = viewLimit # km
    self.interest = interest
    self.powerLimit = powerLimit
  # eventType = [0, 1]
  # 0 - orbit correction
  # 1 - comms object determination

  # Check if the satellite in on its track
  # Return true if on track
  # Return false otherwise, and call direThruster()
  def checkOrbitPath():
    return

  # Given deltaV vector or position offset, determine what and how thruster fires to correct orbit
  # Output the thruster number and duration of the firing to orbit progator
  def fireThruster():
    return
  
  # Given a list of requested objects, determine which objects are in view of current object
  # Return a list of objects in view
  def objectsInView(self, objectsList, stk):
    objectsinView = []
    for obj in objectsList:
      if obj.distance < self.viewLimit:
        objectsinView.append(obj)
    return objectsinView

  # Calculate a link budget of paired satellite between starting and target given transmitting and receiving commos objects
  # Return a list of turples with link budget in form of (start to satellite, satellite to target)
  def calcLinkBudget(start, target, objectsList):
    LB = []
    for obj in objectsList:
      LB.append((LinkBudg.calcPower(start, obj), LinkBudg.calcPower(obj, target)))
    return LB

  # Calculate a DOP given transmitting and receiving commos objects
  # Return a list of DOP of pairs of 2 satellites
  def calcDOP(self, start, target, objectsList):
    DOP = []
    pairs = list(itertools.combinations(objectsList, 2))
    for p in pairs:
      pair = [p[0], p[1]]
      DOP.append((pair, dil_of_pre().calcDOP(pair, self.interest)))
    return DOP

  # Determine what communication obejcts to be created based on optimized link budget and DOP
  # Return a pair of satallite with best DOP considering only power limit
  # Return an empty list if all satallites in view are not able to provide sufficient service with power given
  def finalCommsObjects(self, start, target, objectsList):
    objInView = objectsInView(objectsList, stk)
    DOP = calcDOP(self, start, target, objInView)
    minDOP = DOP.sort(key=lambda x:x[1])
    for pair in minDOP:
      if calcLinkBudget(start,target,pair[0]) < self.powerLimit:
        EventLogger.info("COMMS", "Built communication between " + start + " and " + target + " with satellites #" + str(pair[0][0].name) + " and #" + str(pair[0][1].name))
        return pair[0]
    EventLogger().warning("COMMS", "All available linkages over budget")
    return []
