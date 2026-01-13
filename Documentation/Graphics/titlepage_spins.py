import bpy
import numpy as np
import math as math
import os
import mathutils
import blender_essentials as be
import blender_colors as bc
import vectors as v

r_0 = 0.5

alphaStep = 0.13
tripCol = (1.0,0.0,0.0,1.0)

# ================================================
# =============== Define materials ===============
# ================================================
def createMaterials():
    # surface spin 
    color,alpha = bc.hex_to_rgb(bc.tugreen),0.2
    material = bpy.data.materials.new(name="Sshell")
    material.diffuse_color = (*color,0.1)
    be.render_transparent(material,color,alpha)
    material = bpy.data.materials.new(name="Sspin")
    material.diffuse_color = (*color,1.0)

    for material in bpy.data.materials: # remove existing materials to update them
        material.specular_intensity = 0.5

def addLight():
    # create light datablock, set attributes
    light_data = bpy.data.lights.new(name="LightSource", type='AREA')
    light_data.energy = 700
    light_data.shape = 'SQUARE'
    light_data.size = 10.0

    # create new object with our light datablock
    #light_object = bpy.data.objects.new(name="LightSource1", object_data=light_data)
    #bpy.context.collection.objects.link(light_object)
    #bpy.context.view_layer.objects.active = light_object
    #light_object.location = (-10.07, 5.03, 2.35)
    #euler1 = -76.4
    #euler2 = -58.5
    #euler3 = 52.3
    #light_object.rotation_euler = (euler1*np.pi/180., euler2*np.pi/180., euler3*np.pi/180.)

    # create new object with our light datablock
    light_object = bpy.data.objects.new(name="LightSource2", object_data=light_data)
    bpy.context.collection.objects.link(light_object)
    bpy.context.view_layer.objects.active = light_object
    light_object.location = (0,-4,7.5)#(2.36, -5, 2.63)
    euler1 = 45#44
    euler2 = 0#59.1
    euler3 = 0#-19.3
    light_object.rotation_euler = (euler1*np.pi/180., euler2*np.pi/180., euler3*np.pi/180.)

    # update scene, if needed
    dg = bpy.context.evaluated_depsgraph_get()
    dg.update()

def init():
    be.clearScene()
    be.clearMaterials()
    createMaterials()
    addLight()
    scene = bpy.data.scenes["Scene"]
    cam_rotation = [71.6,0,0] #[59.4, 0, -25.1]
    cam_location = [0,-3.62,2.02] #[-1.59,-3.73,2.6]
    be.positionCamera(scene,cam_rotation,cam_location)

init()

np.random.seed(27) #56

S_size = 0.12

# positions:
S1 = np.array([0,0,0])
S2 = np.array([0.32,0,0.32])
S3 = np.array([0.5,0,0.75])
ensemble = [S1,S2,S3]

sizes = S_size*np.array((1,1.1,1.2))
for S, size in zip(ensemble,sizes):
    SShell = be.SimpleSphere("Sshell",size)
    SSpin = be.Spin("Sspin",size)
    SShell.draw(S)
    SSpin.draw(S,"random")

    

# saving 
bpy.context.scene.render.engine = 'BLENDER_EEVEE' #BLENDER_EEVEE BLENDER_WORKBENCH

res = [900,670]
bpy.context.scene.render.resolution_x = res[0] #1920
bpy.context.scene.render.resolution_y = res[1]  #1080
perc = 400
bpy.context.scene.render.resolution_percentage = perc



#bpy.context.scene.display.shading.light = 'MATCAP'
bpy.context.scene.display.shading.light = 'FLAT'
bpy.context.scene.display.shading.color_type = 'TEXTURE'
bpy.data.scenes["Scene"].render.film_transparent = True
bpy.context.scene.render.image_settings.color_mode = 'RGBA'
filename = r"/home/timo/Schreibtisch/PhD_documents/THESIS/Testenvironment/Graphics/" + "titlepage_spins"
bpy.context.scene.render.filepath = filename
bpy.ops.render.render(write_still = True)

#cropping
scissors_pct = (0.45,0.45,0.75,0.92)
resolution = [perc/100*res[0],perc/100*res[1]]
be.crop_percentage(filename,scissors_pct,resolution=resolution)



# glassy behavior with eevee:
#bpy.context.scene.eevee.use_ssr = True
#bpy.context.scene.eevee.use_ssr_refraction = True
#material.use_screen_refraction = True
#principled_bsdf.inputs["Roughness"].default_value = 0
#principled_bsdf.inputs["Transmission"].default_value = 100