#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PV4FoamReader'
model2_testOpenFOAM = GetActiveSource()
#model2_testOpenFOAM = PV4FoamReader(FileName='/home/caio/OpenFOAM/model2_test/model2_test.OpenFOAM')
model2_testOpenFOAM.MeshParts = ['porosity - cellZone', 'inletchannel - cellZone', 'outletchannel - cellZone']
model2_testOpenFOAM.VolumeFields = ['Con', 'p', 'U']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on model2_testOpenFOAM
#model2_testOpenFOAM.PartArrayStatus = ['internalMesh', '0', 'wall - group', '0', 'back - patch', '0', 'front - patch', '0', 'inlet - patch', '0', 'outlet - patch', '0', 'porosityWall - patch', '0', 'porosity - cellZone', '1', 'inletchannel - cellZone', '1', 'outletchannel - cellZone', '1', 'porosity - faceZone', '0', 'inletchannel - faceZone', '0', 'outletchannel - faceZone', '0']
#model2_testOpenFOAM.MeshParts = ['porosity - cellZone', 'inletchannel - cellZone', 'outletchannel - cellZone']
#model2_testOpenFOAM.VolFieldArrayStatus = ['Con', '1', 'p', '1', 'U', '1']
#model2_testOpenFOAM.VolumeFields = ['Con', 'p', 'U']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1625, 860]

# show data in view
model2_testOpenFOAMDisplay = Show(model2_testOpenFOAM, renderView1)
# trace defaults for the display properties.
model2_testOpenFOAMDisplay.ColorArrayName = [None, '']
model2_testOpenFOAMDisplay.ScalarOpacityUnitDistance = 0.0006636809274881448

# reset view to fit data
renderView1.ResetCamera()

# set scalar coloring
ColorBy(model2_testOpenFOAMDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
model2_testOpenFOAMDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
vtkBlockColorsLUT.InterpretValuesAsCategories = 1
vtkBlockColorsLUT.Annotations = ['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11']
vtkBlockColorsLUT.ActiveAnnotatedValues = ['0', '1', '2']
vtkBlockColorsLUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.63, 0.63, 1.0, 0.67, 0.5, 0.33, 1.0, 0.5, 0.75, 0.53, 0.35, 0.7, 1.0, 0.75, 0.5]

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# create a new 'Slice'
slice1 = Slice(Input=model2_testOpenFOAM)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.0, 9.313225746154785e-10, -5.296897143125534e-09]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1)

# Properties modified on model2_testOpenFOAM
model2_testOpenFOAM.VolumeFields = ['Con', 'p', 'U']
model2_testOpenFOAM.PointFields = []

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.ColorArrayName = [None, '']

# hide data in view
Hide(model2_testOpenFOAM, renderView1)

# set scalar coloring
ColorBy(slice1Display, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'Con'))

# rescale color and/or opacity maps used to include current data range
slice1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Con'
conLUT = GetColorTransferFunction('Con')
conLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 2.5, 0.865003, 0.865003, 0.865003, 5.0, 0.705882, 0.0156863, 0.14902]
conLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Con'
conPWF = GetOpacityTransferFunction('Con')
conPWF.Points = [0.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]
conPWF.ScalarRangeInitialized = 1

# create a new 'Threshold'
threshold1 = Threshold(Input=slice1)
threshold1.Scalars = ['POINTS', 'Con']
threshold1.ThresholdRange = [0.0, 5.0]

# Properties modified on threshold1
threshold1.ThresholdRange = [0.5, 5.0]

# show data in view
threshold1Display = Show(threshold1, renderView1)
# trace defaults for the display properties.
threshold1Display.ColorArrayName = ['POINTS', 'Con']
threshold1Display.LookupTable = conLUT
threshold1Display.ScalarOpacityUnitDistance = 0.00024211868649258222

# hide data in view
Hide(slice1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# current camera placement for renderView1
renderView1.CameraPosition = [0.0, 9.313225746154785e-10, 0.13693568986975505]
renderView1.CameraFocalPoint = [0.0, 9.313225746154785e-10, -5.296897143125534e-09]
renderView1.CameraParallelScale = 0.03544156586348279

# save screenshot
SaveScreenshot('/home/caio/OpenFOAM/model2_test/00.png', magnification=1, quality=100, view=renderView1)

# get animation scene
animationScene1 = GetAnimationScene()

animationScene1.GoToNext()

# current camera placement for renderView1
renderView1.CameraPosition = [0.0, 9.313225746154785e-10, 0.13693568986975505]
renderView1.CameraFocalPoint = [0.0, 9.313225746154785e-10, -5.296897143125534e-09]
renderView1.CameraParallelScale = 0.03544156586348279

# save screenshot
SaveScreenshot('/home/caio/OpenFOAM/model2_test/01.png', magnification=1, quality=100, view=renderView1)

for i in range(2,123):
    animationScene1.GoToNext()

    # current camera placement for renderView1
    renderView1.CameraPosition = [0.0, 9.313225746154785e-10, 0.1131699907086005]
    renderView1.CameraFocalPoint = [0.0, 9.313225746154785e-10, -5.296897143125534e-09]
    renderView1.CameraParallelScale = 0.03544156586348279

    # save screenshot
    SaveScreenshot('/home/caio/OpenFOAM/model2_test/0'+str(i)+'.png', magnification=1, quality=100, view=renderView1)
