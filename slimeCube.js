
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								MASSSPRING DEFINITION  									///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Vec3 {
	constructor( x, y, z ) { this.init(x,y,z); }
	init( x, y, z ) { this.x=x; this.y=y; this.z=z; }
	copy ( ) { return new Vec3( this.x, this.y, this.z ); }
	set  (v) { this.x =v.x; this.y =v.y; this.z =v.z; }
	inc  (v) { this.x+=v.x; this.y+=v.y; this.z+=v.z; }
	dec  (v) { this.x-=v.x; this.y-=v.y; this.z-=v.z; }
	scale(f) { this.x*=f; this.y*=f; this.z*=f; }
	add  (v) { return new Vec3( this.x+v.x, this.y+v.y, this.z+v.z ); }
	sub  (v) { return new Vec3( this.x-v.x, this.y-v.y, this.z-v.z ); }
	dot  (v) { return this.x*v.x + this.y*v.y + this.z*v.z; }
	cross(v) { return new Vec3( this.y*v.z-this.z*v.y, this.z*v.x-this.x*v.z, this.x*v.y-this.y*v.x ); }
	mul  (f) { return new Vec3( this.x*f, this.y*f, this.z*f ); }
	div  (f) { return new Vec3( this.x/f, this.y/f, this.z/f ); }
	len2 ( ) { return this.dot(this); }
	len  ( ) { return Math.sqrt(this.len2()); }
	unit ( ) { return this.div(this.len()); }
	normalize() {
		var l = this.len();
		this.x /= l;
		this.y /= l;
		this.z /= l;
	}
	trans(m) {
		return {
			x: m[0]*this.x + m[4]*this.y + m[ 8]*this.z + m[12],
			y: m[1]*this.x + m[5]*this.y + m[ 9]*this.z + m[13],
			z: m[2]*this.x + m[6]*this.y + m[10]*this.z + m[14],
			w: m[3]*this.x + m[7]*this.y + m[11]*this.z + m[15]
		};
	}
}	
function ToVec3(a) { return new Vec3(a[0],a[1],a[2]); }

class SlimeCube {

	constructor(gl)
	{
		this.gl = gl;

		this.gravity = new Vec3( 0, -2.0, 0 );
		this.mass = .115;
		this.structStiffness = 4.85;
		this.shearStiffness = 4;
		this.bendStiffness = 5;
		this.slider_structDamping = 1.7;
		this.slider_shearDamping  = 1.7;
		this.slider_bendDamping   = 1.15;
		
		this.forwardForceVal=0.2; 
		this.leftForceVal=0.2; 
		this.rightForceVal=0.2; 
		this.backwardForceVal=0.2;
		this.jumpForceVal=0.5;

		this.grad = 0; //serve per il gradiente di colore quando i lcubo mangia
		this.stop_grad_change = true;
		this.restitution = .8;

		this.u_EATSKULL = false;
		this.u_grad = 0;
		this.setMesh(document.getElementById("box.obj").text);

		


		

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								MOVEMENT EVENT LISTENER 								///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
		// Add event listener for key pressure WASD
		this.applyForwardForce=false;
		this.applyLeftForce = false;
		this.applyRightForce = false;
		this.applyBackwardForce = false;
		this.applyJumpForce = false;
		document.addEventListener('keydown', (event) => {
			switch (event.keyCode) {
				case  87: // W (Forward)
				this.applyForwardForce = true;
				break;
				case 65: // A (Left)
				this.applyLeftForce = true;
				break;
				case  68: // D (Right)
				this.applyRightForce = true;
				break;
				case  83: // S (Backward)
				this.applyBackwardForce = true;
				break;
				case 32: // Space (Jump )
				this.applyJumpForce = true; 
				break;
			}
		});
		// Add event listener for  key release
		document.addEventListener('keyup', (event) => {
			switch (event.keyCode) {
				case  87: // W (Forward)
				this.applyForwardForce = false;
				break;
				case 65: // A (Left)
				this.applyLeftForce = false;
				break;
				case  68: // D (Right)
				this.applyRightForce = false;
				break;
				case  83: // S (Backward)
				this.applyBackwardForce = false;
				break;
				case 32: // Space (Jump )
				this.applyJumpForce = false; 
				break;
			}
		});
	}
	
	setMesh( objdef )
		{
			this.obj = parseOBJ( objdef);
			//questo è il parser degl ihomework , che mi da la face e un modo per fare l'update delle normals 
			
			this.mesh = new ObjMesh;
			this.mesh.parse( objdef); 
			this.mesh.computeNormals();
			var box = this.mesh.getBoundingBox();
			
			var size = [
				(box.max[0]-box.min[0])/2,
				(box.max[1]-box.min[1])/2,
				(box.max[2]-box.min[2])/2
			];
	
		
			var maxSize = Math.max( size[0], size[1], size[2] );
			// var scale = 0.4/maxSize;
			this.scale = 1/maxSize;
			// console.log("scale",this.scale, maxSize);
	
			this.parts = this.obj.geometries.map(({data}) => {
			
				const bufferInfo = webglUtils.createBufferInfoFromArrays(this.gl, data);
				return {
				material: {
					u_diffuse: [Math.random(), Math.random(), Math.random(), 1],
				},
				bufferInfo,
				};
			});
	
	
			this.extents = getGeometriesExtents(this.obj.geometries);
			this.range = m4.subtractVectors(this.extents.max, this.extents.min);
			// amount to move the object so its center is at the origin
			this.objOffset = m4.scaleVector(
				m4.addVectors(
					this.extents.min,
					m4.scaleVector(this.range, 0.5)),
				-1);
			this.reset();
			this.initSprings();
			
			
		}
		reset()
		{
			this.pos = Array( this.mesh.vpos.length );
			for ( var i=0; i<this.pos.length; ++i ) this.pos[i] = ToVec3( this.mesh.vpos[i] );
			this.vel = Array( this.pos.length );
			for ( var i=0; i<this.vel.length; ++i ) this.vel[i] = new Vec3(0,0,0);
			this.nrm = Array( this.mesh.norm.length );
			for ( var i=0; i<this.nrm.length; ++i ) this.nrm[i] = ToVec3( this.mesh.norm[i] );
			this.buffers = this.mesh.getVertexBuffers();
			
		}
	
	
	cubeTexture()
	{
		// var furImage = document.getElementById("fur");
		var img = document.getElementById('cube-texture');
		meshDrawer.setTexture( img );
		meshDrawer.showTexture(true);
		meshDrawer.isSlime(true);
					
	}
	updateGradient()
	{
		this.grad+=0.035;
		if (this.grad>=0.7){
			this.grad=0;
			this.stop_grad_change = false;
		}
		return this.grad;
	}


	

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								SPRINGS HANDLING	  									///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// l'obietttivo è creare 3 liste di molle che siano
	// - Structural : Connect every node to its 6 direct neighbours
					// • Node (i,j,k) connected to
					// – (i+1,j,k), (i-1,j,k), (i,j-1,k), (i,j+1,k), (i,j,k-1), (i,j,k+1)
					// (for surface nodes, some of these neighbors might not exists)
					// • Structural springs establish the basic structure
					// of the jello cube

	// - Shear :       • Disallow excessive shearing
					// • Prevent the cube from distorting
					// • Every node (i,j,k)
					// connected to its diagonal
					// neighbors
	// - Bend :    Prevent the cube from folding over
					// • Every node connected
					// to its second neighbor
					// in every direction
					// (6 connections per node,
					// unless surface node)
	identifyNeighbors(pos, threshold) {
	const neighbors = {}; // Dictionary to store neighbors for each vertex (index as key)
	for (let i = 0; i < pos.length; i++) {
		const vertex1 = pos[i];
		const neighborIndices = {
			immediateNeighbors:[],
			diagonalNighboors:[],
			twoDistNeighboors:[]
		};
		for (let j = 0; j < pos.length; j++) {
			if (i !== j) { // Skip comparing a vertex to itself
				const vertex2 = pos[j];
				const distance = this.euclideanDistance(vertex1, vertex2);
				// console.log("[DEBUG 0] distance ", distance);
				// console.log("[DEBUG 0] threshold ", threshold);
				if (distance <= threshold) {
					// console.log("[DEBUG 1] ", j , " in immediate neighbor of " ,i );
				neighborIndices.immediateNeighbors.push(j); // Immediate neighbor
				} else if (distance <=0.0001+(Math.pow(3, 0.5) * threshold)) { // Adjust range for diagonal neighbors (optional)
					// console.log("[DEBUG 2] ", j , " in diagonal  neighbor of " ,i);
					// console.log("[DEBUG 2] , distance ",distance ,"<= (Math.pow(3, 0.5) * threshold) " , (Math.pow(3, 0.5) * threshold));

				neighborIndices.diagonalNighboors.push(j); // Diagonal neighbor
				}
				else if (distance <= 0.001+ ( 2 * (Math.pow(3, 0.5) * threshold)))
				{ // Adjust range for 2-dist  neighbors (optional)
					// console.log("[DEBUG 3] ", j , " in 2-dist neighbor of " ,i);
					// console.log("[DEBUG 3] , distance ",distance ,"<= 2 * (Math.pow(3, 0.5) * threshold) " ,2 * (Math.pow(3, 0.5) * threshold));
				neighborIndices.twoDistNeighboors.push(j); //2-dsit neighbor
				}

			}
		}
		// console.log("neighborIndices",neighborIndices);
		neighbors[i] = neighborIndices;
	}
	return neighbors;
	}

	// Assuming you have a separate function to calculate euclidean distance
	euclideanDistance(p1, p2) {

		var subVector = p2.sub(p1);
		var dist = subVector.len();
		return dist;
	}
	calculateRestLength(p1,p2)
	{
		return p1.sub(p2).len();
	}
	minCubeSpringLen(pos) {
		// Initialize minimum distance with a large value
		let minDistance = Infinity;

		// Loop through all vertices (excluding the diagonal itself)
		for (let i = 0; i < pos.length - 1; i++) {
			const vertex1 = pos[i];
			// Loop through remaining vertices (avoiding self-comparison and already compared pairs)
			for (let j = i + 1; j < pos.length; j++) {
			const vertex2 = pos[j];
			const distance = this.euclideanDistance(vertex1, vertex2); // Use your existing distance calculation logic
			minDistance = Math.min(minDistance, distance); // Update minimum if current distance is smaller
			}
		}

		return minDistance;
	}
	getCubeFacesIndices(vertices) {
	let minXPos = Infinity, maxXPos = -Infinity;
	let minXNeg = Infinity, maxXNeg = -Infinity;
	let minYPos = Infinity, maxYPos = -Infinity;
	let minYNeg = Infinity, maxYNeg = -Infinity;
	let minZPos = Infinity, maxZPos = -Infinity;
	let minZNeg = Infinity, maxZNeg = -Infinity;

	for (let vertex of vertices) {
		if (vertex.x >= 0) {
		minXPos = Math.min(minXPos, vertex.x);
		maxXPos = Math.max(maxXPos, vertex.x);
		} else {
		minXNeg = Math.min(minXNeg, vertex.x);
		maxXNeg = Math.max(maxXNeg, vertex.x);
		}

		if (vertex.y >= 0) {
		minYPos = Math.min(minYPos, vertex.y);
		maxYPos = Math.max(maxYPos, vertex.y);
		} else {
		minYNeg = Math.min(minYNeg, vertex.y);
		maxYNeg = Math.max(maxYNeg, vertex.y);
		}

		if (vertex.z >= 0) {
		minZPos = Math.min(minZPos, vertex.z);
		maxZPos = Math.max(maxZPos, vertex.z);
		} else {
		minZNeg = Math.min(minZNeg, vertex.z);
		maxZNeg = Math.max(maxZNeg, vertex.z);
		}
	}

	const sideLengthX = Math.abs(maxXPos - minXPos);
	const sideLengthY = Math.abs(maxYPos - minYPos);
	const sideLengthZ = Math.abs(maxZPos - minZPos);


	// Assuming all sides are equal, use the minimum side length
	const sideLength = Math.min(sideLengthX, sideLengthY, sideLengthZ);
	const cubeVolume = Math.pow(sideLength, 3);
	const faces = [[], [], [], [], [], []]; // Array to store vertex indices for each face
	const threshold = 0.2; // 20% threshold (adjust as needed)

	for (let i = 0; i < vertices.length; i++) {
		const vertex = vertices[i];

		// Check for front/back faces (based on z-coordinate)
		if (vertex.z >= cubeVolume * threshold) {
		faces[0].push(i); // Front face (positive z)
		} else if (vertex.z <= -cubeVolume * threshold) {
		faces[1].push(i); // Back face (negative z)
		}

		// Check for left/right faces (based on x-coordinate)
		if (vertex.x >= cubeVolume * threshold) {
		faces[2].push(i); // Right face (positive x)
		} else if (vertex.x <= -cubeVolume * threshold) {
		faces[3].push(i); // Left face (negative x)
		}

		// Check for top/bottom faces (based on y-coordinate)
		if (vertex.y >= 1) {
		faces[4].push(i); // Top face (positive y)
		} else if (vertex.y <= 0.3 ){ //-cubeVolume * threshold) {
		faces[5].push(i); // Bottom face (negative y)
		}
	}

	return faces;
	}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								INIT SPRINGS 											///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	initSprings() {
		this.structSprings = [];
		this.shearSprings = [];
		this.bendSprings = [];
		var threshold = this.minCubeSpringLen(this.pos);
		this.cubefaces = this.getCubeFacesIndices(this.pos);
		// Get neighbor information for each node
		const neighbors = this.identifyNeighbors(this.pos, threshold);

		// Loop for Structural Springs
		for (const node in neighbors) {
			// console.log("[DEBUG 0 node]", parseInt(node, 10));
			const immediateNeighbors = neighbors[node].immediateNeighbors;
			for (const neighbor of immediateNeighbors) {
				// console.log("[DEBUG 0 neighbor]", neighbor);
				this.structSprings.push({
					p0: parseInt(node, 10), // Parse node to integer (base 10)
					p1: parseInt(neighbor, 10),
					rest: this.calculateRestLength(this.pos[parseInt(node, 10)], this.pos[parseInt(neighbor, 10)]),
					});
				}
		}

		// Loop for Shear Springs
		for (const node in neighbors) {
			const diagonalNighboors = neighbors[node].diagonalNighboors;
			for (const neighbor of diagonalNighboors) {
				this.shearSprings.push({
					p0: parseInt(node, 10), // Parse node to integer (base 10)
					p1: parseInt(neighbor, 10),
					rest: this.calculateRestLength(this.pos[parseInt(node, 10)], this.pos[parseInt(neighbor, 10)]),
					});
				}
		}

		// Loop for Bend Springs
		for (const node in neighbors) {
			const twoDistNeighboors = neighbors[node].twoDistNeighboors;
			for (const neighbor of twoDistNeighboors) {
				this.bendSprings.push({
					p0: parseInt(node, 10), // Parse node to integer (base 10)
					p1: parseInt(neighbor, 10),
					rest: this.calculateRestLength(this.pos[parseInt(node, 10)], this.pos[parseInt(neighbor, 10)]),
					});
				}
		}
	}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								UPDATE MESH												///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	updateMesh()
	{
		function updateBuffer( buffer, faces, verts )
		{
			function addTriangleToBuffer( buffer, bi, vals, i, j, k )
			{
				buffer[bi++] = vals[i].x;
				buffer[bi++] = vals[i].y;
				buffer[bi++] = vals[i].z;
				buffer[bi++] = vals[j].x;
				buffer[bi++] = vals[j].y;
				buffer[bi++] = vals[j].z;
				buffer[bi++] = vals[k].x;
				buffer[bi++] = vals[k].y;
				buffer[bi++] = vals[k].z;
			}
			for ( var i=0, bi=0; i<faces.length; ++i ) {
				var f = faces[i];
				if ( f.length < 3 ) continue;
				addTriangleToBuffer( buffer, bi, verts, f[0], f[1], f[2] );
				bi += 9;
				for ( var j=3; j<f.length; ++j, bi+=9 ) {
					addTriangleToBuffer( buffer, bi, verts, f[0], f[j-1], f[j] );
				}
			}
		}

		// update the position buffer
		updateBuffer( this.buffers.positionBuffer, this.mesh.face, this.pos );
		// console.log("this.parts[0].bufferInfo.attribs.a_position",this.parts[0].bufferInfo.attribs.a_position);
		// updateBuffer( this.parts[0].bufferInfo.attribs.a_position, this.mesh.face, this.pos );

		// update normals
		for ( var i=0; i<this.nrm.length; ++i ) this.nrm[i].init(0,0,0);
		for ( var i=0; i<this.mesh.face.length; ++i ) {
			var f = this.mesh.face[i];
			var nf = this.mesh.nfac[i];
			var v0 = this.pos[ f[0] ];
			for ( var j=1; j<f.length-1; ++j ) {
				var v1 = this.pos[ f[j] ];
				var v2 = this.pos[ f[j+1] ];
				var e0 = v1.sub(v0);
				var e1 = v2.sub(v0);
				var n  = e0.cross(e1);
				n = n.unit();
				this.nrm[ nf[0  ] ].inc(n);
				this.nrm[ nf[j  ] ].inc(n);
				this.nrm[ nf[j+1] ].inc(n);
			}
		}
		for ( var i=0; i<this.nrm.length; ++i ) this.nrm[i].normalize();

		// console.log("this.parts[0].bufferInfo.attribs.a_normal",this.parts[0].bufferInfo.attribs.a_normal);
		// updateBuffer( this.parts[0].bufferInfo.attribs.a_normal, this.mesh.nfac, this.nrm );
		updateBuffer( this.buffers.normalBuffer, this.mesh.nfac, this.nrm );
		
		const mesh_data = {position : this.buffers.positionBuffer,texcoord:  this.buffers.texCoordBuffer,normal:  this.buffers.normalBuffer};
		const bufferInfo_mesh = webglUtils.createBufferInfoFromArrays(this.gl, mesh_data);
		// console.log("bufferInfo_mesh",bufferInfo_mesh);
		this.parts[0].bufferInfo = bufferInfo_mesh;
		
		// // Update the mesh drawer and redraw scene
		// meshDrawer.setMesh( this.buffers.positionBuffer, this.buffers.texCoordBuffer, this.buffers.normalBuffer );
		
		// // pointDrawer.updatePoint();
		// UpdateCubeViewMatrices(this.pos);
		// DrawScene();
		
		
	}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								SIM TIME STEP FUNCTION CALL								///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	simTimeStep()
	{
		
		// Update positions and velocities
		// var timestep = document.getElementById('timestep').value;
		var timestep = 20;
		const dt = timestep / 1000;	// time step in seconds
		const structDamping = this.slider_structDamping * this.structStiffness * dt;
		const shearDamping = this.slider_shearDamping * this.shearStiffness * dt;
		const bendDamping = this.slider_bendDamping * this.bendStiffness * dt;
		
		const skullCentersArray = Array(skullArray.length);
		for ( var i=0; i<skullCentersArray.length; ++i ) skullCentersArray[i] =average_center( skullArray[i].pos);
		


		var eat_return= SimTimeStep( dt, this.pos, this.vel, this.structSprings,this.shearSprings ,this.bendSprings ,
					this.structStiffness, this.shearStiffness, this.bendStiffness,
					structDamping, shearDamping, bendDamping , this.mass, this.gravity, this.restitution ,
					this.applyForwardForce, this.applyLeftForce ,this.applyRightForce,	this.applyBackwardForce, this.applyJumpForce,
					this.cubefaces , this.forwardForceVal, this.leftForceVal, this.rightForceVal, this.backwardForceVal, this.jumpForceVal,
					skullCentersArray);
		if (typeof eat_return.return_skull_i === 'number' ) {
	
			skullArray[eat_return.return_skull_i].stopDrawing();
			
		}
		
		this.ShowEatSkull(eat_return.EAT_SKULL);
		
		

		this.updateMesh();
	}
	

	


 startSimulation()
{
	this.simTimeStep();
}

ShowEatSkull(eat)
{
	if (!eat && !this.stop_grad_change)
	{
		this.stop_grad_change = true;
		this.u_EATSKULL = false;
		
	}
	var grad=0;
	if (eat && this.stop_grad_change){
		var new_grad= this.updateGradient();
		if (this.stop_grad_change){
			grad = new_grad;
		}
		
	}
	if(this.stop_grad_change){
		this.eatSkull( eat, grad );
	}
	else if(eat && !this.stop_grad_change){
		this.eatSkull( eat, 0.9 );
	}
	
}
updateGradient()
{
	this.grad+=0.035;
	if (this.grad>=0.9){
		this.grad=0;
		this.stop_grad_change = false;
	}
	return this.grad;
}


eatSkull( eat_,grad_ )
{
  this.u_EATSKULL = eat_;
  this.u_grad = grad_;
	
}
	
}

function average_center(vertices){
	let centerX = 0;
	let centerY = 0;
	let centerZ = 0;
	let numVertices = 0; // To keep track of the number of vertices
	for (let i = 0; i < vertices.length; i ++) {
		// Assuming your vertex data is an array where each vertex has x, y, and z components stored consecutively
		centerX += vertices[i].x;
		centerY += vertices[i].y;
		centerZ += vertices[i].z;
		numVertices++;
	}
	centerX /= numVertices;
	centerY /= numVertices;
	centerZ /= numVertices;
	return new Vec3(centerX,centerY,centerZ);

}	








///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								skull				  									///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
const obj_numbers = [0,1,2]; // Your 3 numbers

function pickRandomNumber(numbers) {
	const maxIndex = numbers.length - 1; // Last index of the array
	const randomIndex = Math.floor(Math.random() * (maxIndex + 1)); // Random index between 0 and maxIndex
	return numbers[randomIndex];
}
const NUM_SKULL_TILES =7;
const skull_spacing=[];
var skullArray= [];
	class Skull{
		constructor(gl, i,skullBoolean,bone1Boolean,bone2Boolean){
			this.i = i;
			this.gl = gl;
			var rand1  =Math.random();
			var rand2  =Math.random();
			this.spacing = [(1+Math.random()) * 3,(1+Math.random()) * 3];
			// console.log("spcaing for ", this.i, this.spacing)
			if(rand1>0.5){
				this.spacing[0]=-this.spacing[0]
			}
			if(rand2>0.5){
				this.spacing[1]=-this.spacing[1]
			}
			skull_spacing.push(this.spacing);
			this.isSkull = skullBoolean;
			this.isBone1 = bone1Boolean;
			this.isBone2 = bone2Boolean; 
			if (this.isSkull){
                          
                    this.setMesh( document.getElementById('skull.obj').text );                  
			}
			else if (this.isBone1)  {           
                    this.setMesh( document.getElementById('bone1.obj').text );
                   }
			else if (this.isBone2)   {  
                    this.setMesh( document.getElementById('bone2.obj').text );
                    }
                   

			// this.setMesh( document.getElementById('skull.obj').text );
			this.stop_drawing=false;
			this.r;
			this.c;
			if(NUM_SKULL_TILES == 9){
				if(this.i < 3){
					this.r=0;
					this.c=this.i;
				}
				else if(this.i < 6){
					this.r=1;
					this.c=this.i-3;
				}
				else if(this.i < 9){
					this.r=2;
					this.c=this.i-6;
				}
			}
			else if(NUM_SKULL_TILES == 4){
				if(this.i < 2){
					this.r=0;
					this.c=this.i;
				}
				else if(this.i < 4){
					this.r=1;
					this.c=this.i-2;
				}
			}
		}
		stopDrawing()
		{
			this.stop_drawing=true;
		}
		setMesh( objdef )
		{
		this.obj = parseOBJ( objdef);
		//questo è il parser degl ihomework , che mi da la face e un modo per fare l'update delle normals 
		
		this.mesh = new ObjMesh;
		this.mesh.parse( objdef); 
		this.mesh.computeNormals();
		

        this.parts = this.obj.geometries.map(({data}) => {
           
            const bufferInfo = webglUtils.createBufferInfoFromArrays(this.gl, data);
            return {
            material: {
                u_diffuse: [Math.random(), Math.random(), Math.random(), 1],
            },
            bufferInfo,
            };
        });


        this.extents = getGeometriesExtents(this.obj.geometries);
        this.range = m4.subtractVectors(this.extents.max, this.extents.min);
        // amount to move the object so its center is at the origin
        this.objOffset = m4.scaleVector(
            m4.addVectors(
                this.extents.min,
                m4.scaleVector(this.range, 0.5)),
            -1);

	
			this.reset();
		}
		reset()
		{
			this.pos = Array( this.mesh.vpos.length );
			for ( var i=0; i<this.pos.length; ++i ) this.pos[i] = ToVec3( this.mesh.vpos[i] );
			this.vel = Array( this.pos.length );
			for ( var i=0; i<this.vel.length; ++i ) this.vel[i] = new Vec3(0,0,0);
			this.nrm = Array( this.mesh.norm.length );
			for ( var i=0; i<this.nrm.length; ++i ) this.nrm[i] = ToVec3( this.mesh.norm[i] );
			this.buffers = this.mesh.getVertexBuffers();
			
			}
	}
	



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								FLOOR				  									///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	var NUM_FLOOR_TILES = 4;

	var floorArray= [];
	class Floor{
		constructor(gl, i){
			this.gl = gl;
			this.i = i;
			this.setMesh( document.getElementById('floor.obj').text );
			this.r;
			this.c;
			
			if(this.i < 2){
				this.r=0;
				this.c=this.i;
			}
			else if(this.i < 4){
				this.r=1;
				this.c=this.i-2;
			}
			// else {
			// 	//higher step
			// 	this.r = 0;
			// 	this.c = this.i;
			// }
			
		}
		setMesh( objdef )
		{
			this.obj = parseOBJ( objdef);
			//questo è il parser degl ihomework , che mi da la face e un modo per fare l'update delle normals 
			
			this.mesh = new ObjMesh;
			this.mesh.parse( objdef); 
			this.mesh.computeNormals();
			var box = this.mesh.getBoundingBox();
			
			var size = [
				(box.max[0]-box.min[0])/2,
				(box.max[1]-box.min[1])/2,
				(box.max[2]-box.min[2])/2
			];
	
		
			var maxSize = Math.max( size[0], size[1], size[2] );
			// var scale = 0.4/maxSize;
			this.scale = 1/maxSize;
			// console.log("scale",this.scale, maxSize);

			this.parts = this.obj.geometries.map(({data}) => {
			
				const bufferInfo = webglUtils.createBufferInfoFromArrays(this.gl, data);
				return {
				material: {
					u_diffuse: [Math.random(), Math.random(), Math.random(), 1],
				},
				bufferInfo,
				};
			});


			this.extents = getGeometriesExtents(this.obj.geometries);
			this.range = m4.subtractVectors(this.extents.max, this.extents.min);
			// amount to move the object so its center is at the origin
			this.objOffset = m4.scaleVector(
				m4.addVectors(
					this.extents.min,
					m4.scaleVector(this.range, 0.5)),
				-1);
			this.reset();
			
			
			
		}
		reset()
		{
			this.pos = Array( this.mesh.vpos.length );
			for ( var i=0; i<this.pos.length; ++i ) this.pos[i] = ToVec3( this.mesh.vpos[i] );
			this.vel = Array( this.pos.length );
			for ( var i=0; i<this.vel.length; ++i ) this.vel[i] = new Vec3(0,0,0);
			this.nrm = Array( this.mesh.norm.length );
			for ( var i=0; i<this.nrm.length; ++i ) this.nrm[i] = ToVec3( this.mesh.norm[i] );
			this.buffers = this.mesh.getVertexBuffers();
			
		}
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								WALL				  									///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

var NUM_WALL_TILES = 10;

var wallArray= [];
class Wall{
	constructor(gl, i){
		this.gl = gl;
		this.i = i;
		this.setMesh( document.getElementById('wall.obj').text );
		this.r;
		this.c;
		this.rot = 0;
		
		if(this.i < 2){
			this.r=0;
			this.c=this.i;
			this.rot = 0;
		}
		else if(this.i < 4){
			this.r=1;
			this.c=this.i-2;
			this.rot = 0;
		}
		else if(this.i < 6){
			this.r=0;
			this.c=this.i-4;
			this.rot = Math.PI / 2;
		}
		else if(this.i < 8){
			this.r=1;
			this.c=this.i-6;
			this.rot = Math.PI / 2;
		}
		else{
			//corrispettivi di 0 e 5 
			//muri dello step
			if (this.i ==8){ //0
				this.r=0;
				this.c=0;
				this.rot = 0;
			}
			else if( this.i==9){ //5 
				this.r=0;
				this.c=1;
				this.rot = Math.PI / 2;
			}
		}

		
	}
	setMesh( objdef )
	{
		this.obj = parseOBJ( objdef);
		//questo è il parser degl ihomework , che mi da la face e un modo per fare l'update delle normals 
		
		this.mesh = new ObjMesh;
		this.mesh.parse( objdef); 
		this.mesh.computeNormals();
		var box = this.mesh.getBoundingBox();
		
		var size = [
			(box.max[0]-box.min[0])/2,
			(box.max[1]-box.min[1])/2,
			(box.max[2]-box.min[2])/2
		];

	
		var maxSize = Math.max( size[0], size[1], size[2] );
		// var scale = 0.4/maxSize;
		this.scale = 1/maxSize;
		// console.log("scale",this.scale, maxSize);

		this.parts = this.obj.geometries.map(({data}) => {
		
			const bufferInfo = webglUtils.createBufferInfoFromArrays(this.gl, data);
			return {
			material: {
				u_diffuse: [Math.random(), Math.random(), Math.random(), 1],
			},
			bufferInfo,
			};
		});


		this.extents = getGeometriesExtents(this.obj.geometries);
		this.range = m4.subtractVectors(this.extents.max, this.extents.min);
		// amount to move the object so its center is at the origin
		this.objOffset = m4.scaleVector(
			m4.addVectors(
				this.extents.min,
				m4.scaleVector(this.range, 0.5)),
			-1);
		this.reset();
		
		
		
	}
	reset()
	{
		this.pos = Array( this.mesh.vpos.length );
		for ( var i=0; i<this.pos.length; ++i ) this.pos[i] = ToVec3( this.mesh.vpos[i] );
		this.vel = Array( this.pos.length );
		for ( var i=0; i<this.vel.length; ++i ) this.vel[i] = new Vec3(0,0,0);
		this.nrm = Array( this.mesh.norm.length );
		for ( var i=0; i<this.nrm.length; ++i ) this.nrm[i] = ToVec3( this.mesh.norm[i] );
		this.buffers = this.mesh.getVertexBuffers();
		
	}
}




















/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


DEBUGVEL = false;
STOP = false;
DEBUGFORCE = false;
DEBUGSPRING = false;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								AUX FUNCTIONS											///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


accumulated_jump = 0;


accumulated_forward = 0;
forward_deceleration  = 0;
stop_forward_dec=true; //da il segnale di smettere di decelerare 
forward_decelerating = false; // dice che si sta decelerando

accumulated_backward = 0;
backward_deceleration = 0;
stop_backward_dec = true;
backward_decelerating = false;


accumulated_left = 0;
left_deceleration = 0;
stop_left_dec = true;
left_decelerating = false;


accumulated_right= 0;
right_deceleration = 0;
stop_right_dec = true;
right_decelerating = false;

EAT_SKULL = false;
eaten_skull_enter=-1;
eaten_skull_array_exit=[];
////////////////////////////////////////////////////////DAMPING FORCE///////////////////////////////////////////////
function damping_force( p1Pos, p0Pos, spring, tempVelocities, damping){
	// Calculate direction vector (B to A)
	//springvector
	var direction = p1Pos.sub(p0Pos);
	if (DEBUGSPRING){
		console.log("[DAMPING] direction ",direction);
	}
	// Calculate relative velocity
	var relativeVelocity = tempVelocities[spring.p0].sub(tempVelocities[spring.p1]);
	if (DEBUGSPRING){
		
		console.log("[DAMPING] relativeVelocity ",relativeVelocity);
	
	}
	// Project relative velocity onto direction vector
	// var projectedRelVel = direction.cross(relativeVelocity.cross(direction).div(direction.len2()));
	var projectedRelVel = relativeVelocity.cross(direction.unit());
	if (DEBUGSPRING){
		
		console.log("[DAMPING] projectedRelVel ",projectedRelVel);
	}
	// Calculate damping force contribution for this spring
	var springDampingForce = projectedRelVel.mul(-damping).cross(direction.unit());
	if (DEBUGSPRING){// 
		console.log("[DAMPING] projectedRelVel.mul(-damping) ",projectedRelVel.mul(-damping));
		console.log("[DAMPING] direction.unit() ",direction.unit());
		console.log("[DAMPING] springDampingForce ",springDampingForce);
	}

	return springDampingForce;
}

function spring_force(p1Pos, p0Pos, stiffness, spring){
	var springVector = p1Pos.sub(p0Pos);
	var springLength = springVector.len();

	var springForce = springVector.unit().mul(-stiffness * (springLength - spring.rest));
	return springForce;
}
function springforce_and_damping(spring, stiffness,  tempVelocities, damping, tempPosition){
	
		if (DEBUGSPRING){
			console.log( "Spring : ",spring );
		}
		var p0Pos = tempPosition[spring.p0];
		var p1Pos = tempPosition[spring.p1];
		var springForce = spring_force(p1Pos, p0Pos, stiffness, spring);
		var springDampingForce = damping_force(p1Pos, p0Pos, spring, tempVelocities, damping);
		
		
		if (DEBUGSPRING){
			console.log( "springForce : ",springForce );
			console.log( "springDampingForce : ",springDampingForce );
		}
		let p0 = spring.p0,
			p1 = spring.p1;
		return {springForce, springDampingForce, p0, p1} ;
}
function initializeForces(positions) {
	return Array(positions.length).fill(new Vec3(0, 0, 0));
  }

function avarage_center(vertices){
	let centerX = 0;
	let centerY = 0;
	let centerZ = 0;
	let numVertices = 0; // To keep track of the number of vertices
	for (let i = 0; i < vertices.length; i ++) {
		// Assuming your vertex data is an array where each vertex has x, y, and z components stored consecutively
		centerX += vertices[i].x;
		centerY += vertices[i].y;
		centerZ += vertices[i].z;
		numVertices++;
	}
	centerX /= numVertices;
	centerY /= numVertices;
	centerZ /= numVertices;
	return new Vec3(centerX,centerY,centerZ);

}

function euclideanDistance(p1, p2) {

	var subVector = p2.sub(p1);
	var dist = subVector.len();
	return dist;
}

function computeCubeSize(positions) {
	// Minimum expected positions (8 vertices)
	if (positions.length < 8) {
	  return null; // Not enough points for a cube
	}
  
	// Create a set of unique positions
	const uniquePositions = new Set(positions);
  
	// Minimum and maximum coordinates to define a bounding box
	let minX = Infinity, maxX = -Infinity;
	let minY = Infinity, maxY = -Infinity;
	let minZ = Infinity, maxZ = -Infinity;
  
	for (const pos of uniquePositions) {
	  minX = Math.min(minX, pos.x);
	  maxX = Math.max(maxX, pos.x);
	  minY = Math.min(minY, pos.y);
	  maxY = Math.max(maxY, pos.y);
	  minZ = Math.min(minZ, pos.z);
	  maxZ = Math.max(maxZ, pos.z);
	}
  
	// Expected size difference between adjacent vertices (tolerance)
	const tolerance = 1e-6;
  
	// Potential cube vertices based on bounding box extremes
	const candidateVertices = [
	  new Vec3(minX, minY, minZ),
	  new Vec3(maxX, minY, minZ),
	  new Vec3(maxX, maxY, minZ),
	  new Vec3(minX, maxY, minZ),
	  new Vec3(minX, minY, maxZ),
	  new Vec3(maxX, minY, maxZ),
	  new Vec3(maxX, maxY, maxZ),
	  new Vec3(minX, maxY, maxZ),
	];
  
	// Verify candidate vertices (check distances and potential duplicates)
	const identifiedVertices = [];
	for (const candidate of candidateVertices) {
	  let foundMatch = false;
	  for (const pos of uniquePositions) {
		if (Math.abs(euclideanDistance(candidate,pos)) <= tolerance) {
		  if (foundMatch) {
			throw new Error("Duplicate vertex found during cube identification.");
		  }
		  foundMatch = true;
		  identifiedVertices.push(pos);
		  break; // Skip further checks for this candidate if a close match is found
		}
	  }
	  if (!foundMatch) {
		// Candidate vertex not found in actual positions (not a perfect cube?)
		return null;
	  }
	}
  
	// Calculate pairwise distances and store in a set
	const allDistances = new Set();
	for (let i = 0; i < identifiedVertices.length; i++) {
	  for (let j = i + 1; j < identifiedVertices.length; j++) {
		allDistances.add(euclideanDistance(identifiedVertices[i],identifiedVertices[j]));
	  }
	}
  
	// Calculate the average edge size
	const averageSize = Array.from(allDistances).reduce((sum, distance) => sum + distance, 0) / allDistances.size;
  
	return averageSize;
  }












///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								SIMULATION FUNCTION 									///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function SimTimeStep( dt, positions, velocities, structSprings,shearSprings,bendSprings, 
	structStiffness, shearStiffness, bendStiffness, 
	structDamping, shearDamping, bendDamping ,
	 particleMass, gravity, restitution ,
	  applyForwardForce  ,applyLeftForce ,applyRightForce, applyBackwardForce, applyJumptForce,
	faceArray, forwardForceVal, leftForceVal, rightForceVal, backwardForceVal, jumpForceVal,
	skullCentersArray)
{

	
	// [TO-DO] Compute the total force of each particle
	var forces = Array(positions.length).fill(new Vec3(0, 0, 0)); // Initialize forces array with zeros
	// const forces = initializeForces(positions);
	const forwardForce = new Vec3( 0, 0, forwardForceVal );
    const leftForce = new Vec3(-leftForceVal, 0, 0);
    const rightForce = new Vec3(rightForceVal, 0, 0);
    const backwardForce = new Vec3(0, 0,-backwardForceVal);
	const jumpForce = new Vec3(0, -jumpForceVal, 0);
	
	


	// Create a copy of the velocity array
	const tempVelocities = [];
	for (var i = 0; i < velocities.length; i++) {
	  tempVelocities.push(velocities[i].copy()); // Use the copy method of your Vec3 class
	}
	const tempPosition = [];
	for (var i = 0; i < positions.length; i++) {
		tempPosition.push(positions[i].copy()); // Use the copy method of your Vec3 class
		tempPosition[i].z +=3;
	}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////								DECELERATION PT.1/3	 									///////////////


////////////////////////////////////////////FORWARD////////////////////////////////////////////////////////////////

//check per movimenti inintenzionali 
	for (var i = 0; i < positions.length; i++) {
		if(!applyForwardForce && faceArray[5].includes(i) && accumulated_forward <= 0 && forward_decelerating && stop_forward_dec){
			//se ci sono movimenti in avanti . li stoppo
			if (tempVelocities[i].z>0){
				tempVelocities[i].z=0;
			}
		}
	}
	//check per la decelerazione
	if(!applyForwardForce){ //viene rilasciato W
		if(stop_forward_dec && !forward_decelerating){// se è true , allora non c'è decelerazione
			stop_forward_dec = false;// si pò calcolare la dec
		}
	}

////////////////////////////////////////////BACKWARD////////////////////////////////////////////////////////////////

	//check per movimenti inintenzionali 
	for (var i = 0; i < positions.length; i++) {
		if(!applyBackwardForce && faceArray[5].includes(i) && accumulated_backward <= 0 && backward_decelerating && stop_backward_dec){
			//se ci sono movimenti in avanti . li stoppo
			if (tempVelocities[i].z<0){
				tempVelocities[i].z=0;
			}
		}
	}
	//check per la decelerazione
	if(!applyBackwardForce){ //viene rilasciato W
		if(stop_backward_dec && !backward_decelerating){// se è true , allora non c'è decelerazione
			stop_backward_dec = false;// si pò calcolare la dec
		}
	}



////////////////////////////////////////////LEFT////////////////////////////////////////////////////////////////

	//check per movimenti inintenzionali 
	for (var i = 0; i < positions.length; i++) {
		if(!applyLeftForce && faceArray[5].includes(i) && accumulated_left <= 0 && left_decelerating && stop_left_dec){
			//se ci sono movimenti in avanti . li stoppo
			if (tempVelocities[i].x<0){
				tempVelocities[i].x=0;
			}
		}
	}
	//check per la decelerazione
	if(!applyLeftForce){ //viene rilasciato W
		if(stop_left_dec && !left_decelerating){// se è true , allora non c'è decelerazione
			stop_left_dec = false;// si pò calcolare la dec
		}
	}



////////////////////////////////////////////RIGHT////////////////////////////////////////////////////////////////

	//check per movimenti inintenzionali 
	for (var i = 0; i < positions.length; i++) {
		if(!applyRightForce && faceArray[5].includes(i) && accumulated_right <= 0 && right_decelerating && stop_right_dec){
			//se ci sono movimenti in avanti . li stoppo
			if (tempVelocities[i].x>0){
				tempVelocities[i].x=0;
			}
		}
	}
	//check per la decelerazione
	if(!applyRightForce){ //viene rilasciato W
		if(stop_right_dec && !right_decelerating){// se è true , allora non c'è decelerazione
			stop_right_dec = false;// si pò calcolare la dec
		}
	}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								SPRINGS FORCCES		 									///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (DEBUGSPRING){
	console.log("structSprings", structSprings);
	console.log("shearSprings", shearSprings);
	console.log("bendSprings", bendSprings);
	}
	///structSprings
	for (var i = 0; i < structSprings.length; i++) {
		var spring = structSprings[i];
		var spring_and_dam = springforce_and_damping(
			spring,
			structStiffness,  
			tempVelocities, 
			structDamping, 
			tempPosition);
		const  structSpringForce = spring_and_dam.springForce, 
			   structSpringDampingForce = spring_and_dam.springDampingForce, 
			   p0 = spring_and_dam.p0, 
			   p1 = spring_and_dam.p1;
		if(DEBUGSPRING){
			console.log("p0,p1: ", p0, p1 ,"structSpringForce : ",structSpringForce," structSpringDampingForce : ", structSpringDampingForce );
		}
		////////////////ADD FORCES///////////////////
		forces[p0] = forces[p0].sub(structSpringForce).sub(structSpringDampingForce);
		forces[p1] = forces[p1].add(structSpringForce).add(structSpringDampingForce);

	}
	//shear springs
	for (var i = 0; i < shearSprings.length; i++) {
		var spring = shearSprings[i];
		var spring_and_dam = springforce_and_damping(
			spring,
			shearStiffness,  
			tempVelocities, 
			shearDamping, 
			tempPosition);
		
		const  shearSpringForce = spring_and_dam.springForce, 
			   shearSpringDampingForce = spring_and_dam.springDampingForce, 
			   p0 = spring_and_dam.p0, 
			   p1 = spring_and_dam.p1;
		if(DEBUGSPRING){
			console.log("p0,p1: ", p0, p1 ,"shearSpringForce : ",shearSpringForce," shearSpringDampingForce : ", shearSpringDampingForce );
		}
		////////////////ADD FORCES///////////////////
		forces[p0] = forces[p0].sub(shearSpringForce).sub(shearSpringDampingForce);
		forces[p1] = forces[p1].add(shearSpringForce).add(shearSpringDampingForce);

	}
	//bend springs
	for (var i = 0; i < bendSprings.length; i++) {
		var spring = bendSprings[i];
		var spring_and_dam = springforce_and_damping(
			spring,
			bendStiffness,  
			tempVelocities, 
			bendDamping, 
			tempPosition);
		const  bendSpringForce = spring_and_dam.springForce, 
			   bendSpringDampingForce = spring_and_dam.springDampingForce, 
			   p0 = spring_and_dam.p0, 
			   p1 = spring_and_dam.p1;
		
		
		if(DEBUGSPRING){
			console.log("p0,p1: ", p0, p1 ,"bendSpringForce : ",bendSpringForce," bendSpringDampingForce : ", bendSpringDampingForce );
		}
		////////////////ADD FORCES///////////////////
		forces[p0] = forces[p0].sub(bendSpringForce).sub(bendSpringDampingForce);
		forces[p1] = forces[p1].add(bendSpringForce).add(bendSpringDampingForce);

	}
		
	
	
	if (DEBUGFORCE) {
		// console.log(" ---beforeup all forces : ",tempForces );
		console.log(" ---beforeup all forces : ",forces );
	}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								CONTROL MOVEMENT FORCES									///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	var cube_center = avarage_center(tempPosition);
	// var cubeSize  = computeCubeSize(positions);
	// console.log("cubeSize",cubeSize)


	for (var i = 0; i < positions.length; i++) {
		if (DEBUGFORCE) {
			// console.log(i, " ---beforeup force : ",tempForces[i] );
			console.log(i, " ---beforeup force : ",forces[i] );
		}

		
		/////////////////////////////////////////////////////////////////////////////
		//////////////				DECELERATION PT.2/3							/////
			
		/////////////////////////////////////////FORWARD//////////////////////////////////////////////////////
		if (applyForwardForce){
			if(faceArray[5].includes(i)){
				//face 0 is forward
				forces[i]  = forces[i].add(forwardForce);
				accumulated_forward+=backwardForceVal;
				forward_decelerating = false;

				// console.log("forwardForce",forwardForce,"\naccumulated_forward",accumulated_forward, "\nforces[i]", forces[i] );

			}
			
		}
		else if (accumulated_forward>0){
			if(faceArray[5].includes(i)){
				if (!forward_decelerating && !stop_forward_dec){ 
					forward_deceleration = accumulated_forward*0.1/100
					forward_decelerating = true;
				}
				forces[i]  = forces[i].add(new Vec3(0, 0,-forward_deceleration));
				accumulated_forward-=forward_deceleration;
			
				

			}
			
		}

// 		accumulated_backward = 0;
// backward_deceleration = 0;
// stop_backward_dec = true;
// backward_decelerating = false;
		/////////////////////////////////////////BACKWARD//////////////////////////////////////////////////////
		if (applyBackwardForce){
			//face 1 is back 
			if(faceArray[5].includes(i)){
				forces[i]  = forces[i].add(backwardForce);
				accumulated_backward+=forwardForceVal;
				backward_decelerating = false;
			}
		}
		else if (accumulated_backward>0){
			if(faceArray[5].includes(i)){
				if (!backward_decelerating && !stop_backward_dec){ 
					backward_deceleration = accumulated_backward*0.1/100
					backward_decelerating = true;
				}
				forces[i]  = forces[i].add(new Vec3(0, 0,backward_deceleration));
				accumulated_backward-=backward_deceleration;
			}
			
		}

// accumulated_left = 0;
// left_deceleration = 0;
// stop_left_dec = true;
// left_decelerating = false;

		/////////////////////////////////////////LEFT//////////////////////////////////////////////////////
		if (applyLeftForce){
			//face 1 is back 
			if(faceArray[5].includes(i)){
				forces[i]  = forces[i].add(leftForce);
				accumulated_left+=leftForceVal;
				left_decelerating = false;
			}
		}
		else if (accumulated_left>0){
			if(faceArray[5].includes(i)){
				if (!left_decelerating && !stop_left_dec){ 
					left_deceleration = accumulated_left *0.1/100
					left_decelerating = true;
				}
				forces[i]  = forces[i].add(new Vec3( left_deceleration,0, 0,));
				accumulated_left-=left_deceleration;
			}
			
		}

// accumulated_right= 0;
// right_deceleration = 0;
// stop_right_dec = true;
// right_decelerating = false;
		/////////////////////////////////////////RIGHT//////////////////////////////////////////////////////
		if (applyRightForce){
			//face 1 is back 
			if(faceArray[5].includes(i)){
				forces[i]  = forces[i].add(rightForce);
				accumulated_right+=rightForceVal;
				right_decelerating = false;
			}
		}
		else if (accumulated_right>0){
			if(faceArray[5].includes(i)){
				if (!right_decelerating && !stop_right_dec){ 
					right_deceleration = accumulated_right *0.1/100
					right_decelerating = true;
				}
				forces[i]  = forces[i].add(new Vec3(- right_deceleration,0, 0,));
				accumulated_right-=right_deceleration;
			}
			
		}
		

		/////////////////////////////////JUMP////////////////////////////////////////////
	
		
	
		if (applyJumptForce){
			if(faceArray[4].includes(i)){
				//face4 is up
				forces[i]  = forces[i].add(jumpForce);
				accumulated_jump+=jumpForceVal/100;
			}
		}
		else if (accumulated_jump>0){
			if(faceArray[4].includes(i)){
				forces[i]  = forces[i].add(new Vec3(0,accumulated_jump,0));
			}
		}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								GRAVITY FORCE											///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		forces[i]  = forces[i].add(gravity.mul(particleMass));

		
		
		if (DEBUGFORCE) {
			console.log(i, " ---G forces: ",gravity.mul(particleMass) )
			console.log(i, " ---damping forces: ",tempVelocities[i].mul(-damping));
			console.log(i, " ---afterup forces : ",forces[i] );
		}

	}
	if (!applyJumptForce){
		accumulated_jump = 0; //reset the jump
	}



	if (DEBUGFORCE) {
		console.log(" AFTER ---FORCES : ",forces );
	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								UPDATE VELOCITIES AND POSITION							///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Update velocities using Semi-implicit Euler
	for (var i = 0; i < positions.length; i++) {
		if(DEBUGVEL){
			console.log(i, "updatating velocities :  tempVelocities[i]->",tempVelocities[i]);
			console.log(i, "forces[i]: ",forces[i] , "vel to add ", forces[i].mul(dt).div(particleMass) );
		}

		tempVelocities[i] = tempVelocities[i].add(forces[i].mul(dt).div(particleMass));
		

		/////////////////////////////////////////////////////////////////////////////
		//////////////				DECELERATION PT.3/3			
		
		//////////////////////////////////FORWARD///////////////////////////////////////////
		if(faceArray[5].includes(i)){											
			if(accumulated_forward <= 0 && forward_decelerating && !stop_forward_dec){  
				forward_deceleration = 0;												
				stop_forward_dec = true;										
				tempVelocities[i].z =0;											
			}																	
		}	
		
		

// 		accumulated_backward = 0;
// backward_deceleration = 0;
// stop_backward_dec = true;
// backward_decelerating = false;
		//////////////////////////////////BACKWARD///////////////////////////////////////////
		if(faceArray[5].includes(i)){											
			if(accumulated_backward <= 0 && backward_decelerating && !stop_backward_dec){  
				backward_deceleration = 0;												
				stop_backward_dec = true;										
				tempVelocities[i].z =0;											
			}																	
		}	

// accumulated_left = 0;
// left_deceleration = 0;
// stop_left_dec = true;
// left_decelerating = false;
		//////////////////////////////////LEFT///////////////////////////////////////////
		if(faceArray[5].includes(i)){											
			if(accumulated_left <= 0 && left_decelerating && !stop_left_dec){  
				left_deceleration = 0;												
				stop_left_dec = true;										
				tempVelocities[i].x =0;											
			}																	
		}	


// accumulated_right= 0;
// right_deceleration = 0;
// stop_right_dec = true;
// right_decelerating = false;
		//////////////////////////////////RIGHT///////////////////////////////////////////
		if(faceArray[5].includes(i)){											
			if(accumulated_right <= 0 && right_decelerating && !stop_right_dec){  
				right_deceleration = 0;												
				stop_right_dec = true;										
				tempVelocities[i].x =0;											
			}																	
		}	



		/////////////////////////////////////////////////////////////////////////////


		if(DEBUGVEL){
			console.log(i, " final  tempVelocities  : ",tempVelocities[i] );
			STOP = true;
		}
	}


	// Update positions
	for (var i = 0; i < positions.length; i++) {
		// positions[i] = positions[i].add(velocities[i].mul(dt));
		tempPosition[i] = tempPosition[i].add(tempVelocities[i].mul(dt));
	}
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								EATING HANDLING  										///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	function min_dist_skull_aux(array){
		//array [ [k, dist], [], ...]
		var min_dist=100;
		var min_k = -1
		for(var i = 0; i< array.length;i++){
			var k = array[i][0];
			var dist = array[i][1];
			min_dist = Math.min(min_dist,dist);
			if (min_dist == dist){
				min_k = k;
			}
			
		}
		return {min_dist, min_k}
	}


	var epsilon = 0.0001; // Define a small epsilon value
	var bounce = false; // true to make it bouunce , false to five 0 velocity when collide 
	var return_skull_i;
	if(skullCentersArray.length==eaten_skull_array_exit.length){
		//allora tutti gli skull sono gia stati mangiati 
		EAT_SKULL = false;
	}
	else{
		var min_dist_array = [];
		for(var k =0; k<skullCentersArray.length;k++){
			
			if(!(eaten_skull_array_exit.includes(k))){
				var up_spacing = 0;
				if( skull_spacing[k][0] < 0 && skull_spacing[k][1]<0){
					up_spacing  =4;
				}
				var skullCenter = skullCentersArray[k].add(new Vec3(skull_spacing[k][0],up_spacing,skull_spacing[k][1])); 
				var dist_from_skull= euclideanDistance(cube_center,skullCenter);
				min_dist_array.push([k,dist_from_skull]);
				
			}
			
		}
		if(min_dist_array.length>0){
			var {min_dist, min_k}  =  min_dist_skull_aux(min_dist_array);
			// console.log("min_dist",min_dist, "\n min_k",min_k)
			if (min_dist> 1.2){
				EAT_SKULL = false;

			}
			else if (min_dist< 1){
				EAT_SKULL = true;
				eaten_skull_enter = min_k;
			}
		}
		if (eaten_skull_enter != -1  && !eaten_skull_array_exit.includes(k)&&!EAT_SKULL)
			{
				eaten_skull_array_exit.push(eaten_skull_enter);
				return_skull_i = eaten_skull_enter;
				eaten_skull_enter =-1;
				EAT_SKULL = true;
			}
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////																						///////////////
////////////								WALLS COLLISION   									///////////////
////////////																						///////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for (var i = 0; i < positions.length; i++) {
		var pos = tempPosition[i];
		const groundDepth = 1; // Adjust this value to define the ground plane's z-coordinate
		
		// Check for bottom collision (y-axis)

		if (pos.y + epsilon < -groundDepth) {
			// DEBUGSPRING = true;
			// DEBUGFORCE = true;
			// DEBUGVEL = true;
	
			var penetrationDepth = Math.abs(pos.y + groundDepth);
			tempPosition[i].y = -groundDepth + restitution * penetrationDepth;
			if (bounce){tempVelocities[i].y = tempVelocities[i].y * -restitution;}
			else{tempVelocities[i].y = 0.0;}

			
		}
  
	  // Check for top collision (y-axis)
	  if (pos.y + epsilon > groundDepth*25) {
			var penetrationDepth =  Math.abs(pos.y + groundDepth); // Optional boundary correction
			tempPosition[i].y = -groundDepth + restitution * penetrationDepth;
			if (bounce){tempVelocities[i].y = tempVelocities[i].y * -restitution;}
			else{tempVelocities[i].y = 0.0;}
	  }
	  
	//[TEST]   ADDED  *2 TO ACCOUNT FOR LARGER BOX . IS FOR TESTING , TO BE REMOVED 
	  // Check for left collision (x-axis)

	  const mult_size = 10 ;
	  if (pos.x - epsilon < -groundDepth*mult_size) {
			if (bounce){tempVelocities[i].x = tempVelocities[i].y * -restitution;}
			else{tempVelocities[i].x = 0.0;}
		}
  
	  // Check for right collision (x-axis)
	  if (pos.x + epsilon > groundDepth*mult_size) {
			if (bounce){tempVelocities[i].x = tempVelocities[i].y * -restitution;}
			else{tempVelocities[i].x = 0.0;}
	  }
	  
	  // Check for back collision (z-axis)
	  if (pos.z - epsilon < -groundDepth*mult_size) {
			if (bounce){tempVelocities[i].z = tempVelocities[i].y * -restitution;}
			else{tempVelocities[i].z = 0.0;}
	  }
  
	  // Check for front collision (z-axis)
	  if (pos.z + epsilon >groundDepth*mult_size) {
			if (bounce){tempVelocities[i].z = tempVelocities[i].y * -restitution;}
			else{tempVelocities[i].z = 0.0;}
	  }

	  //check for left top square 
	  if (pos.x + epsilon <0 && pos.z + epsilon < 0) {
		// console.log(pos.y );

		// 0 is  the groud depth of elevated step
		if (pos.y + epsilon <3 && pos.x  < -0.5 && pos.z  < -0.5){
			var penetrationDepth = Math.abs(pos.y + 0);
			tempPosition[i].y = 3 ;//+ restitution * penetrationDepth;
			tempVelocities[i].y = 0.0;
		}
		if (pos.y <0 ){	
			console.log("[left top square  collision (z-axis)]",pos);
		tempVelocities[i].x = 0.0;
		tempVelocities[i].z = 0.0;
	}
	}
	

}
	
	// Replace the original velocity array with the updated copy
	velocities.length = 0; // Clear the original array (optional)
	for (var i = 0; i < tempVelocities.length; i++) {
	  velocities.push(tempVelocities[i]);
	  if(DEBUGVEL){
		console.log(i,"replace velocities");
		console.log("original before replace:" , velocities[i]) ;
		console.log("temp :" , tempVelocities[i]) ;
		console.log("original after replace:" , velocities[i]) ;
		
	  }
	}
	
	positions.length = 0; // Clear the original array (optional)
	for (var i = 0; i < tempPosition.length; i++) {
		positions.push(new Vec3(tempPosition[i].x,tempPosition[i].y,tempPosition[i].z-3));
	  if(DEBUGSPRING){
		console.log(i,"replace positions");
		// console.log("original before replace:" , velocities[i]) ;
		// console.log("temp :" , tempVelocities[i]) ;
		console.log("original after replace:" , positions[i]) ;
		
	  }
	}
	
	if(STOP){
			DEBUGSPRING = false;
			STOP = false;
			DEBUGFORCE = false;
			DEBUGVEL = false;
		}
	return {EAT_SKULL, return_skull_i}




	


}




