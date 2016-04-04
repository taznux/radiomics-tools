import os
import vtk
import qt
import ctk
import slicer
from EditorLib import *

#########################################################
#
#
comment = """

  InterpolateROIsEffect is a subclass of MorphologyEffect
  to Interpolate a layer of pixels from a labelmap

# TODO :
"""
#
#########################################################

#
# InterpolateROIsEffectOptions - see Effect for superclasses
#

class InterpolateROIsEffectOptions(EditorLib.MorphologyEffectOptions):
  """ InterpolateROIsEffect-specfic gui
  """

  def __init__(self, parent=0):
    super(InterpolateROIsEffectOptions,self).__init__(parent)

  def __del__(self):
    super(InterpolateROIsEffectOptions,self).__del__()

  def create(self):
    super(InterpolateROIsEffectOptions,self).create()

    self.minimumSigma = 0.01
    self.maximumSigma = 1.5

    self.sigmaFrame = qt.QFrame(self.frame)
    self.sigmaFrame.setLayout(qt.QHBoxLayout())
    self.frame.layout().addWidget(self.sigmaFrame)
    self.widgets.append(self.sigmaFrame)
    self.sigmaLabel = qt.QLabel("Sigma:", self.sigmaFrame)
    self.sigmaLabel.setToolTip("Set the sigma of the gaussian filter")
    self.sigmaFrame.layout().addWidget(self.sigmaLabel)
    self.widgets.append(self.sigmaLabel)
    self.sigmaSpinBox = slicer.qMRMLSpinBox(self.sigmaFrame)
    self.sigmaSpinBox.objectName = 'SpinBox_Sigma'
    self.sigmaSpinBox.setToolTip("Set the sigma of the gaussian filter")
    #self.sigmaSpinBox.quantity = "length"
    # QFlags not wrapped in python. Equivalent to Prefix | Suffix
    # See qMRMLSpinBox for more details.
    self.sigmaSpinBox.unitAwareProperties = 0x01 | 0x02
    self.sigmaSpinBox.decimals = 2
    self.sigmaSpinBox.minimum = self.minimumSigma
    self.sigmaSpinBox.maximum = self.maximumSigma
    self.sigmaSpinBox.singleStep = self.minimumSigma
    self.sigmaSpinBox.setMRMLScene(slicer.mrmlScene)
    self.sigmaFrame.layout().addWidget(self.sigmaSpinBox)
    self.widgets.append(self.sigmaSpinBox)

    self.sigma = ctk.ctkDoubleSlider(self.frame)
    self.sigma.objectName = 'DoubleSlider_Sigma'
    self.sigma.minimum = self.minimumSigma
    self.sigma.maximum = self.maximumSigma
    self.sigma.orientation = 1
    self.sigma.singleStep = self.minimumSigma
    self.frame.layout().addWidget(self.sigma)
    self.widgets.append(self.sigma)


    self.apply = qt.QPushButton("Apply", self.frame)
    self.apply.objectName = self.__class__.__name__ + 'Apply'
    self.apply.setToolTip("Interpolate current label")
    self.frame.layout().addWidget(self.apply)
    self.widgets.append(self.apply)

    EditorLib.HelpButton(self.frame, "Use this tool to remove pixels from the boundary of the current label.")

    self.connections.append( (self.sigma, 'valueChanged(double)', self.onSigmaValueChanged) )
    self.connections.append( (self.sigmaSpinBox, 'valueChanged(double)', self.onSigmaSpinBoxChanged) )
    self.connections.append( (self.apply, 'clicked()', self.onApply) )

    # Add vertical spacer
    self.frame.layout().addStretch(1)

    self.parameterNode.SetParameter("MorphologyEffect,neighborMode","8")
    self.parameterNode.SetParameter("InterpolateROIsEffect,sigma","1.0")

  def destroy(self):
    super(InterpolateROIsEffectOptions,self).destroy()

  def onApply(self):
    logic = InterpolateROIsEffectLogic(self.editUtil.getSliceLogic())
    logic.undoRedo = self.undoRedo
    fill = int(self.parameterNode.GetParameter('MorphologyEffect,fill'))
    neighborMode = self.parameterNode.GetParameter('MorphologyEffect,neighborMode')
    iterations = int(self.parameterNode.GetParameter('MorphologyEffect,iterations'))
    sigma = float(self.parameterNode.GetParameter('InterpolateROIsEffect,sigma'))
    logic.close(fill,neighborMode,iterations)
    logic.smooth(fill,sigma,iterations)

  # note: this method needs to be implemented exactly as-is
  # in each leaf subclass so that "self" in the observer
  # is of the correct type
  def updateParameterNode(self, caller, event):
    node = self.editUtil.getParameterNode()
    if node != self.parameterNode:
      if self.parameterNode:
        node.RemoveObserver(self.parameterNodeTag)
      self.parameterNode = node
      self.parameterNodeTag = node.AddObserver(vtk.vtkCommand.ModifiedEvent, self.updateGUIFromMRML)

  def setMRMLDefaults(self):
    super(InterpolateROIsEffectOptions,self).setMRMLDefaults()
    disableState = self.parameterNode.GetDisableModifiedEvent()
    self.parameterNode.SetDisableModifiedEvent(1)
    if self.parameterNode.GetParameter("InterpolateROIsEffect,sigma") == '':
      self.parameterNode.SetParameter("InterpolateROIsEffect,sigma","1.0")
    self.parameterNode.SetDisableModifiedEvent(disableState)

  def updateGUIFromMRML(self,caller,event):
    if self.parameterNode.GetParameter("InterpolateROIsEffect,sigma") == '':
      # don't update if the parameter node has not got all values yet
      return
    super(InterpolateROIsEffectOptions,self).updateGUIFromMRML(caller,event)
    self.disconnectWidgets()
    sigma = float(self.parameterNode.GetParameter("InterpolateROIsEffect,sigma"))
    self.sigmaSpinBox.setValue(sigma)
    self.sigma.setValue(sigma)
    self.connectWidgets()

  def onSigmaValueChanged(self,value):
    self.sigmaSpinBox.setValue(self.sigma.value)
    self.updateMRMLFromGUI()

  def onSigmaSpinBoxChanged(self,value):
    self.sigma.setValue(self.sigmaSpinBox.value)
    self.updateMRMLFromGUI()

  def updateMRMLFromGUI(self):
    disableState = self.parameterNode.GetDisableModifiedEvent()
    self.parameterNode.SetDisableModifiedEvent(1)
    super(InterpolateROIsEffectOptions,self).updateMRMLFromGUI()
    self.parameterNode.SetParameter("InterpolateROIsEffect,sigma",str(self.sigmaSpinBox.value))
    self.parameterNode.SetDisableModifiedEvent(disableState)
    if not disableState:
      self.parameterNode.InvokePendingModifiedEvent()


#
# InterpolateROIsEffectTool
#

class InterpolateROIsEffectTool(EditorLib.MorphologyEffectTool):
  """
  One instance of this will be created per-view when the effect
  is selected.  It is responsible for implementing feedback and
  label map changes in response to user input.
  This class observes the editor parameter node to configure itself
  and queries the current view for background and label volume
  nodes to operate on.
  """

  def __init__(self, sliceWidget):
    super(InterpolateROIsEffectTool,self).__init__(sliceWidget)

  def cleanup(self):
    """
    call superclass to clean up actors
    """
    super(InterpolateROIsEffectTool,self).cleanup()

#
# InterpolateROIsEffectLogic
#

class InterpolateROIsEffectLogic(EditorLib.MorphologyEffectLogic):
  """
  This class contains helper methods for a given effect
  type.  It can be instanced as needed by an InterpolateROIsEffectTool
  or InterpolateROIsEffectOptions instance in order to compute intermediate
  results (say, for user feedback) or to implement the final
  segmentation editing operation.  This class is split
  from the InterpolateROIsEffectTool so that the operations can be used
  by other code without the need for a view context.
  """

  def __init__(self,sliceLogic):
    super(InterpolateROIsEffectLogic,self).__init__(sliceLogic)

  def close(self,fill,neighborMode,iterations):

    dilater = slicer.vtkImageErode()
    eroder = slicer.vtkImageErode()
    if vtk.VTK_MAJOR_VERSION <= 5:
      dilater.SetInput( self.getScopedLabelInput() )
    else:
      dilater.SetInputData( self.getScopedLabelInput() )

    if vtk.VTK_MAJOR_VERSION <= 5:
      eroder.SetInput( dilater.GetOutput() )
    else:
      eroder.SetInputData( dilater.GetOutput() )
    eroder.SetOutput( self.getScopedLabelOutput() )

    dilater.SetForeground( fill )
    dilater.SetBackground( self.editUtil.getLabel() )

    eroder.SetForeground( self.editUtil.getLabel() )
    eroder.SetBackground( fill )

    if neighborMode == '8':
      dilater.SetNeighborTo8()
      eroder.SetNeighborTo8()
    elif neighborMode == '4':
      dilater.SetNeighborTo4()
      eroder.SetNeighborTo4()
    else:
      # TODO: error feedback from effect logic?
      # bad neighbor mode - silently use default
      print('Bad neighborMode: %s' % neighborMode)

    for i in xrange(iterations):
      # TODO: $this setProgressFilter dilater "Interpolate ($i)"
      dilater.Update()
      eroder.Update()

    self.applyScopedLabel()
    dilater.SetOutput( None )
    eroder.SetOutput( None )

  def smooth(self,fill,sigma,iterations):
    volumeNode = self.editUtil.getLabelVolume()
    if not volumeNode:
      return

    parameters={}

    parameters["labelToSmooth"] = self.editUtil.getLabel()
    parameters["gaussianSigma"] = sigma
    parameters["numberOfIterations"] = iterations*5
    parameters["inputVolume"] = volumeNode.GetID()
    parameters["outputVolume"] = volumeNode.GetID()

    print parameters

    smoother = slicer.modules.labelmapsmoothing

    slicer.cli.run(smoother, None, parameters, wait_for_completion=True)
    #self.applyScopedLabel()

#
# The InterpolateROIsEffect class definition
#

class InterpolateROIsEffectExtension(EditorLib.MorphologyEffect):
  """Organizes the Options, Tool, and Logic classes into a single instance
  that can be managed by the EditBox
  """

  def __init__(self):
    # name is used to define the name of the icon image resource (e.g. InterpolateROIsEffect.png)
    self.name = "InterpolateROIsEffect"
    # tool tip is displayed on mouse hover
    self.toolTip = "Interpolate: add boundary pixel layers for labelmap editing"

    self.options = InterpolateROIsEffectOptions
    self.tool = InterpolateROIsEffectTool
    self.logic = InterpolateROIsEffectLogic

""" Test:

sw = slicer.app.layoutManager().sliceWidget('Red')
import EditorLib
pet = EditorLib.InterpolateROIsEffectTool(sw)

"""

#
# InterpolateROIsEffect
#

class InterpolateROIsEffect:
  """
  This class is the 'hook' for slicer to detect and recognize the extension
  as a loadable scripted module
  """
  def __init__(self, parent):
    parent.dependencies = ["Editor"]
    parent.title = "Editor Interpolate Effect"
    parent.categories = ["Developer Tools.Editor Extensions"]
    parent.contributors = ["Wookjin Choi"] # insert your name in the list
    parent.helpText = """ """
    parent.acknowledgementText = """ """


    # TODO:
    # don't show this module - it only appears in the Editor module
    #parent.hidden = True

    # Add this extension to the editor's list for discovery when the module
    # is created.  Since this module may be discovered before the Editor itself,
    # create the list if it doesn't already exist.
    try:
      slicer.modules.editorExtensions
    except AttributeError:
      slicer.modules.editorExtensions = {}
    slicer.modules.editorExtensions['InterpolateROIsEffect'] = InterpolateROIsEffectExtension

#
# InterpolateROIsEffectWidget
#

class InterpolateROIsEffectWidget:
  def __init__(self, parent = None):
    self.parent = parent

  def setup(self):
    # don't display anything for this widget - it will be hidden anyway
    pass

  def enter(self):
    pass

  def exit(self):
    pass


