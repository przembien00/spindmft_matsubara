import bpy
import numpy as np
import math as math
import vectors as v
import os
import mathutils
r5 = lambda value : np.around(value,5) # round values for terminal output
pi = 3.14159265

# ===== MATERIALS =====
def compute_rotation_angles(dir):
    n = v.normalize(dir)
    phi = 0
    if r5(np.abs(n[2])) != 1.: # otherwise the projection does not exist
        n_projection = v.normalize(np.array([n[0],n[1],0])) # normalized projection to XY
        sgn = 2*(n[1]>=0)-1 # sign of ny is required for the angle determination
        phi = sgn*np.arccos(n_projection[0]) + pi/2 # pi/2 has to be added since the cylinder is rotated to the -y axis 
    theta = np.arccos(n[2])
    return phi, theta

# fixed cone cylinder radius ratio and radius lenght ratio
class Spin:
    def __init__(self,sort,size):  
        self.sort=sort
        self.length = size
        self.radius = size * 1/10

    def draw(self,pos,dir):
        if dir == "random":
            rands = np.random.rand(1,3)[0]
            dir = np.array([r-0.5 for r in rands])
            
        # cylinder and cone position
        cyl_pos = pos - self.length/5*v.normalize(dir)
        cone_pos = pos + self.length/2*v.normalize(dir)
        
        # rotation
        phi,theta = compute_rotation_angles(dir) 
        rot = (theta,0,phi)
        
        # draw cylinder
        bpy.ops.mesh.primitive_cylinder_add(radius=self.radius,depth=self.length,
                                            location=cyl_pos,rotation=rot)
        bpy.context.active_object.data.materials.append(bpy.data.materials[self.sort])
        bpy.ops.object.shade_smooth()
        
        # draw cone
        crad = self.radius*4   # cone radius
        clen = self.length*3/4 # cone length
        cone = bpy.ops.mesh.primitive_cone_add(radius1=crad,radius2=0,depth=clen,
                                               location=cone_pos,rotation=rot)
        bpy.context.active_object.data.materials.append(bpy.data.materials[self.sort])
        bpy.ops.object.shade_smooth()

# fixed cone cylinder radius ratio
class Vector:
    def __init__(self,sort,cone_radius):  
        self.sort=sort
        self.cyl = SimpleConnection(sort,cone_radius * 1/4)
        self.cone = SimpleConeNection(sort,cone_radius)

    def draw(self,ri,rj,cyl_length_ratio):
        rk = ri + cyl_length_ratio*(rj-ri)
        self.cyl.draw(ri,rk)
        self.cone.draw(rk,rj)

# fixed cone absolute length
class Vector2:
    def __init__(self,sort,cone_radius):  
        self.sort=sort
        self.cone_radius = cone_radius
        self.cyl = SimpleConnection(sort,cone_radius * 1/5)
        self.cone = SimpleConeNection(sort,cone_radius)

    def draw(self,ri,rj):
        rk = rj - 2.5*self.cone_radius*v.normalize(rj-ri)
        self.cyl.draw(ri,rk)
        self.cone.draw(rk,rj)

class SimpleSphere:
    def __init__(self,sort,radius):
        self.sort = sort
        self.radius = radius

    def draw(self,position):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=self.radius,location=position)
        bpy.context.active_object.data.materials.append(bpy.data.materials[self.sort])
        bpy.ops.object.shade_smooth()

class SimpleAtom:
    def __init__(self,sort):  
        self.sort = sort

    def draw(self, position):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=bpy.data.materials[self.sort]["radius"], location=position)
        bpy.context.active_object.data.materials.append(bpy.data.materials[self.sort])
        bpy.ops.object.shade_smooth()

class SimpleConnection:
    def __init__(self,sort,radius):
        self.sort = sort
        self.radius = radius
    
    def draw(self,ri,rj):
            r_rel = v.dist(ri,rj) # distance vector
            r_center = (ri + rj)/2.0 # position of center of bond
            phi, theta = compute_rotation_angles(r_rel)
            bpy.ops.mesh.primitive_cylinder_add(radius=self.radius,depth=v.norm(r_rel),location=r_center,rotation=(theta,0,phi))
            bpy.context.active_object.data.materials.append(bpy.data.materials[self.sort])
            bpy.ops.object.shade_smooth()

class SimpleConeNection:
    def __init__(self,sort,radius):
        self.sort = sort
        self.radius = radius
    
    def draw(self,ri,rj):
            r_rel = v.dist(ri,rj) # distance vector
            r_center = (ri + rj)/2.0 # position of center of bond
            phi, theta = compute_rotation_angles(r_rel)
            bpy.ops.mesh.primitive_cone_add(radius1=self.radius,radius2=0,depth=v.norm(r_rel),location=r_center,rotation=(theta,0,phi))
            bpy.context.active_object.data.materials.append(bpy.data.materials[self.sort])
            bpy.ops.object.shade_smooth()

# ===== OTHERS =====
def positionCamera(scene,rotation,location,is_ortho=False):
    bpy.ops.object.camera_add()
    cam = bpy.data.objects['Camera']
    cam.rotation_mode = 'XYZ'
    if is_ortho:
        cam.data.type ="ORTHO"
    scene.camera = cam
    scene.camera.rotation_euler = (pi / 180.0)*np.array(rotation)
    scene.camera.location.x = location[0]
    scene.camera.location.y = location[1]
    scene.camera.location.z = location[2]
    area = next(area for area in bpy.context.screen.areas if area.type == 'VIEW_3D')
    area.spaces[0].region_3d.view_perspective = 'CAMERA'
    
def clearScene():
    #objs.remove(objs["LightSource"], do_unlink=True)  
    bpy.ops.object.select_all(action='DESELECT')
    for obj in bpy.context.scene.objects:
        obj.select_set(True)  # Select the object
        bpy.ops.object.delete()  # Delete the object

def clearMaterials():
    for material in bpy.data.materials: # remove existing materials to update them
        bpy.data.materials.remove(material)
    
def render_transparent(material,color,alpha):
    material.use_nodes = True
    principled_bsdf = material.node_tree.nodes.get("Principled BSDF")
    principled_bsdf.inputs["Base Color"].default_value = (*color,1)
    principled_bsdf.inputs["Alpha"].default_value = alpha
    material.blend_method = "HASHED"
    material.shadow_method = "HASHED"
    return principled_bsdf

# scissors : left upper right lower
def crop(filename,scissors):
    from PIL import Image
    image_path = filename + ".png"
    image = Image.open(image_path)
    cropped_image = image.crop(scissors)
    cropped_image.save(image_path)

# scissors : left upper right lower
std_resolution = [1920,1080] # x, y
def crop_percentage(filename,scissors_pct,resolution=std_resolution):
    from PIL import Image
    image_path = filename + ".png"
    image = Image.open(image_path)
    scissors = [np.rint(scissors_pct[0]*resolution[0]), np.rint(scissors_pct[1]*resolution[1]), np.rint(scissors_pct[2]*resolution[0]), np.rint(scissors_pct[3]*resolution[1])]
    cropped_image = image.crop(scissors)
    cropped_image.save(image_path)