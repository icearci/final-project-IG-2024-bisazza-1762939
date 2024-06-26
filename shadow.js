

'use strict';
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////                                                                         /////////
///////////////                             MAIN                                        /////////
///////////////                                                                         /////////
/////////////////////////////////////////////////////////////////////////////////////////////////
function main() {
  // Get A WebGL context
  /** @type {HTMLCanvasElement} */
  const canvas = document.getElementById("canvas");
  const gl = canvas.getContext('webgl');
  if (!gl) {
    return;
  }

  const ext = gl.getExtension('WEBGL_depth_texture');
  if (!ext) {
    return alert('need WEBGL_depth_texture');  // eslint-disable-line
  }

  // setup GLSL programs
  const textureProgramInfo = webglUtils.createProgramInfo(gl, ['vertex-shader-3d', 'fragment-shader-3d']);
  const colorProgramInfo = webglUtils.createProgramInfo(gl, ['color-vertex-shader', 'color-fragment-shader']);

 
  const cubeLinesBufferInfo = webglUtils.createBufferInfoFromArrays(gl, {
    position: [
      -1, -1, -1,
       1, -1, -1,
      -1,  1, -1,
       1,  1, -1,
      -1, -1,  1,
       1, -1,  1,
      -1,  1,  1,
       1,  1,  1,
    ],
    indices: [
      0, 1,
      1, 3,
      3, 2,
      2, 0,

      4, 5,
      5, 7,
      7, 6,
      6, 4,

      0, 4,
      1, 5,
      3, 7,
      2, 6,
    ],
  });

  /////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////                                                                         /////////
  ///////////////                             OBJ initialization                          /////////
  ///////////////                                                                         /////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  const slime  =  new SlimeCube(gl);
  for (var i=0; i<NUM_SKULL_TILES;i++)
	{
			const randomChoice =  pickRandomNumber(obj_numbers);
			// console.log("Randomly chosen number:", randomChoice);
			var skullBoolean;
			var bone1Boolean;
			var bone2Boolean; 
			if (randomChoice == 0){
                          
				skullBoolean=true;
				bone1Boolean=false;
				bone2Boolean=false;
			}
			else if (randomChoice == 1)  {           
				skullBoolean=false;
				bone1Boolean=true;
				bone2Boolean=false;
			}
			else if (randomChoice == 2)   {  
				skullBoolean=false;
				bone1Boolean=false;
				bone2Boolean=true;
			}
		skullArray.push(new Skull(gl,i,skullBoolean
			, bone1Boolean
			, bone2Boolean));

	}

for (var i=0; i<NUM_FLOOR_TILES;i++){
  floorArray.push(new Floor(gl,i));
}


for (var i=0; i<NUM_WALL_TILES;i++){
  wallArray.push(new Wall(gl,i));
}

  /////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////                                                                         /////////
  ///////////////                             TEXTURES                                    /////////
  ///////////////                                                                         /////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  const img = document.getElementById("cube-texture")
  const slime_texture = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, slime_texture);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGB, gl.RGB, gl.UNSIGNED_BYTE, img);
  gl.generateMipmap(gl.TEXTURE_2D);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

 

  // make a 8x8 checkerboard texture
  const checkerboardTexture = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, checkerboardTexture);
  gl.texImage2D(
      gl.TEXTURE_2D,
      0,                // mip level
      gl.LUMINANCE,     // internal format
      8,                // width
      8,                // height
      0,                // border
      gl.LUMINANCE,     // format
      gl.UNSIGNED_BYTE, // type
      new Uint8Array([  // data
        0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC,
        0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF,
        0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC,
        0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF,
        0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC,
        0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF,
        0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC,
        0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF, 0xCC, 0xFF,
      ]));
  gl.generateMipmap(gl.TEXTURE_2D);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

  /////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////                                                                         /////////
  ///////////////                            depthtexuter                                 /////////
  ///////////////                                                                         /////////
  /////////////////////////////////////////////////////////////////////////////////////////////////

  //we create a texture then a framebuffer and a attach the texture to the framebuffer as a DEPTH_ATTACHMENT
  const depthTexture = gl.createTexture();
  const depthTextureSize = 512;
  gl.bindTexture(gl.TEXTURE_2D, depthTexture);
  gl.texImage2D(
      gl.TEXTURE_2D,      // target
      0,                  // mip level
      gl.DEPTH_COMPONENT, // internal format
      depthTextureSize,   // width
      depthTextureSize,   // height
      0,                  // border
      gl.DEPTH_COMPONENT, // format
      gl.UNSIGNED_INT,    // type
      null);              // data
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

  const depthFramebuffer = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, depthFramebuffer);
  gl.framebufferTexture2D(
      gl.FRAMEBUFFER,       // target
      gl.DEPTH_ATTACHMENT,  // attachment point
      gl.TEXTURE_2D,        // texture target
      depthTexture,         // texture
      0);                   // mip level

  // create a color texture of the same size as the depth texture
  // done because of webgl internal reasons (in order to work on all browser)
  const unusedTexture = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, unusedTexture);
  gl.texImage2D(
      gl.TEXTURE_2D,
      0,
      gl.RGBA,
      depthTextureSize,
      depthTextureSize,
      0,
      gl.RGBA,
      gl.UNSIGNED_BYTE,
      null,
  );
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

  // attach it to the framebuffer
  gl.framebufferTexture2D(
      gl.FRAMEBUFFER,        // target
      gl.COLOR_ATTACHMENT0,  // attachment point
      gl.TEXTURE_2D,         // texture target
      unusedTexture,         // texture
      0);                    // mip level

  function degToRad(d) {
    return d * Math.PI / 180;
  }

  const settings = {
    run: 0,
    cameraX: 6,
    cameraY: 9, //5
    cameraZ: -18,//-7,
    posX: 2.5,
    posY: 8.0,//4.8
    posZ: 4.3,
    targetX: 1.4 , // 2.5,
    targetY: 5, //0
    targetZ: 3.5,
    projWidth: 1,
    projHeight: 1,
    perspective: true,
    fieldOfView: 140,
    bias: -0.006,
  
  };
  webglLessonsUI.setupUI(document.querySelector('#ui'), settings, [
    { type: 'slider',   key: 'run',    min: 0, max: 1, change: render, precision: 1, step: 1, },
    { type: 'slider',   key: 'cameraX',    min: -20, max: 20, change: render, precision: 2, step: 0.001, },
    { type: 'slider',   key: 'cameraY',    min:   -20, max: 20, change: render, precision: 2, step: 0.001, },
    { type: 'slider',   key: 'cameraZ',    min:   -20, max: 20, change: render, precision: 2, step: 0.001, },
    { type: 'slider',   key: 'posX',       min: -20, max: 20, change: render, precision: 2, step: 0.001, },
    { type: 'slider',   key: 'posY',       min:   -20, max: 20, change: render, precision: 2, step: 0.001, },
    { type: 'slider',   key: 'posZ',       min:   -20, max: 20, change: render, precision: 2, step: 0.001, },
    { type: 'slider',   key: 'targetX',    min: -10, max: 10, change: render, precision: 2, step: 0.001, },
    { type: 'slider',   key: 'targetY',    min:   0, max: 20, change: render, precision: 2, step: 0.001, },
    { type: 'slider',   key: 'targetZ',    min: -10, max: 20, change: render, precision: 2, step: 0.001, },
    { type: 'slider',   key: 'projWidth',  min:   0, max:  2, change: render, precision: 2, step: 0.001, },
    { type: 'slider',   key: 'projHeight', min:   0, max:  2, change: render, precision: 2, step: 0.001, },
    { type: 'checkbox', key: 'perspective', change: render, },
    { type: 'slider',   key: 'fieldOfView', min:  1, max: 179, change: render, },
    { type: 'slider',   key: 'bias',       min:  -0.01, max: 0.00001, change: render, precision: 4, step: 0.0001, },
  ]);

  const fieldOfViewRadians = degToRad(60);

  
 //start the simulation
  if(settings.run == 1){slime.startSimulation(slime);}


  /////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////                                                                         /////////
  ///////////////                             DRAW SCENE                                  /////////
  ///////////////                                                                         /////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  function drawScene(
      projectionMatrix,
      cameraMatrix,
      textureMatrix,
      lightWorldMatrix,
      time,
      programInfo) {
    // Make a view matrix from the camera matrix.
    const viewMatrix = m4.inverse(cameraMatrix);

    gl.useProgram(programInfo.program);

    // set uniforms that are the same for both the sphere and plane
    // note: any values with no corresponding uniform in the shader
    // are ignored.
    webglUtils.setUniforms(programInfo, {
      u_view: viewMatrix,
      u_projection: projectionMatrix,
      u_bias: settings.bias,
      u_textureMatrix: textureMatrix,
      u_projectedTexture: depthTexture,
      u_shininess: 150,
      u_innerLimit: Math.cos(degToRad(settings.fieldOfView / 2 - 10)),
      u_outerLimit: Math.cos(degToRad(settings.fieldOfView / 2)),
      u_lightDirection: lightWorldMatrix.slice(8, 11).map(v => -v),
      u_lightWorldPosition: [settings.posX, settings.posY, settings.posZ],
      u_viewWorldPosition: cameraMatrix.slice(12, 15),
    });

  
 /////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////                                                                         /////////
  ///////////////                             DRAW SKULLS                                  /////////
  ///////////////                                                                         /////////
  /////////////////////////////////////////////////////////////////////////////////////////////////

    for (var i=0; i<NUM_SKULL_TILES;i++)
    {
      var skull = skullArray[i];
      if(!skullArray[i].stop_drawing){
        for (const {bufferInfo, material} of skull.parts) {
          // calls gl.bindBuffer, gl.enableVertexAttribArray, gl.vertexAttribPointer
          webglUtils.setBuffersAndAttributes(gl, programInfo, bufferInfo);
          // calls gl.uniform
          
          // let u_world = m4.translation(skull.spacing*skull.r, 1,skull.spacing*skull.c);
          var up_spacing = 1;
          if( skull.spacing[0] < 0 && skull.spacing[1]<0){
            up_spacing  =5;
          }

          let u_world = m4.translation(skull.spacing[0], up_spacing,skull.spacing[1]);
          u_world = m4.translate(u_world, ...skull.objOffset);
          u_world = m4.scale(u_world, 0.2, 0.2, 0.2);
        //  u_world = m4.yRotate(u_world, time );  // Rotate around the translated position

          const objCubeUniforms = {
          u_colorMult: [1, 1, 1, 1],  
          u_color: [0, 0, 1, 1],
          u_texture: checkerboardTexture,
          u_world,
          // u_world: m4.translation(6, 1, 0),
          u_diffuse: material.u_diffuse
        };  
          webglUtils.setUniforms(programInfo, objCubeUniforms);
          // calls gl.drawArrays or gl.drawElements
          webglUtils.drawBufferInfo(gl, bufferInfo);
        }
      }

    }
     /////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////                                                                         /////////
  ///////////////                             DRAW FLOOR                                  /////////
  ///////////////                                                                         /////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  for (var i=0; i<NUM_FLOOR_TILES;i++)
  {
    var floor = floorArray[i];
    
      for (const {bufferInfo, material} of floor.parts) {
        // calls gl.bindBuffer, gl.enableVertexAttribArray, gl.vertexAttribPointer
        webglUtils.setBuffersAndAttributes(gl, programInfo, bufferInfo);
        // calls gl.uniform
        
        // let u_world = m4.translation(-15*floor.scale, 0,15*floor.scale);
        // u_world = m4.translate(u_world, ...floor.objOffset);
        // u_world = m4.scale(u_world, 1,1,1);//floor.scale,floor.scale, floor.scale);
        // console.log(u_world)
        var scaling = floor.scale *5;
        var spacing = 100;
        var up_translation = 0;
        if(floor.i == 2){
          up_translation=40;// step sopraelevato
        }
        let u_world = m4.scaling(scaling,scaling,scaling);
        u_world = m4.translate(u_world, 50-spacing*floor.r,up_translation, -50+spacing*floor.c);
        u_world = m4.translate(u_world, ...floor.objOffset);
          

        const objCubeUniforms = {
        u_colorMult: [0.5, 0.5, 0.5, 1],  
        u_color: [0, 0, 1, 1],
        u_texture: checkerboardTexture,
        u_world,
        // u_world: m4.translation(6, 1, 0),
        u_diffuse: material.u_diffuse
      };  
        webglUtils.setUniforms(programInfo, objCubeUniforms);
        // calls gl.drawArrays or gl.drawElements
        webglUtils.drawBufferInfo(gl, bufferInfo);
      }
    

  }


     /////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////                                                                         /////////
  ///////////////                             DRAW WALL                                  /////////
  ///////////////                                                                         /////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  for (var i=0; i<NUM_WALL_TILES;i++)
  {
    var wall = wallArray[i];
    
      for (const {bufferInfo, material} of wall.parts) {
        // calls gl.bindBuffer, gl.enableVertexAttribArray, gl.vertexAttribPointer
        webglUtils.setBuffersAndAttributes(gl, programInfo, bufferInfo);
        
        
        // console.log(u_world)
        var scaling = wall.scale *5.5;
        var spacing =17.5;
        var deepness = -17.5;
        if(wall.i >7){
          deepness = 0;
        }

        let u_world = m4.scaling(scaling,scaling,scaling); 
        u_world = m4.yRotate(u_world,wall.rot);
        u_world = m4.translate(u_world,-8+spacing*wall.c, 3+6.5*wall.r, deepness);
        u_world = m4.translate(u_world, ...wall.objOffset);
       
 
        const objCubeUniforms = {
        u_colorMult: [0.3, 0.0, 0.0, 1],  
        u_color: [0, 0, 1, 1],
        u_texture: checkerboardTexture,
        u_world,
        // u_world: m4.translation(6, 1, 0),
        u_diffuse: material.u_diffuse
      };  
        webglUtils.setUniforms(programInfo, objCubeUniforms);
        // calls gl.drawArrays or gl.drawElements
        webglUtils.drawBufferInfo(gl, bufferInfo);
      }
    

  }
 /////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////                                                                         /////////
  ///////////////                             DRAW SLIME                                  /////////
  ///////////////                                                                         /////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
    function mixColors(color1, color2, weight) {
      const wr = weight * (color2[0] - color1[0]) + color1[0];
      const wg = weight * (color2[1] - color1[1]) + color1[1];
      const wb = weight * (color2[2] - color1[2]) + color1[2];
      const wa = weight * (color2[3] - color1[3]) + color1[3];
      return [wr, wg, wb, wa];
    }
   // compute the world matrix once since all parts
   // are at the same space.
  
   for (const {bufferInfo, material} of slime.parts) {
     // calls gl.bindBuffer, gl.enableVertexAttribArray, gl.vertexAttribPointer
     webglUtils.setBuffersAndAttributes(gl, programInfo, bufferInfo);
     // calls gl.uniform
   
     var scaling = slime.scale*1;
    let u_world = m4.scaling(scaling,scaling,scaling); 
     u_world = m4.translate(u_world, 0,2, 3);
     u_world = m4.translate(u_world, ...slime.objOffset);
    //  u_world = m4.yRotate(u_world, time );  // Rotate around the translated position

    let u_colorMult= [0.3, 1.0, 0.0, 1];
    let reddish_tint = [1.0, 0.0, 0.0, 1]; // Reddish tint
    if (slime.u_EATSKULL){
      u_colorMult = mixColors(u_colorMult, reddish_tint, slime.u_grad);
    }
     const objCubeUniforms = {
      u_colorMult,  
      u_color: [0, 0, 1, 1],
      u_texture: slime_texture,
      u_world,
      // u_world: m4.translation(6, 1, 0),
      u_diffuse: material.u_diffuse
    };  
     webglUtils.setUniforms(programInfo, objCubeUniforms);
     // calls gl.drawArrays or gl.drawElements
     webglUtils.drawBufferInfo(gl, bufferInfo);
   }


  }

 

  /////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////                                                                         /////////
  ///////////////                             RENDER                                      /////////
  ///////////////                                                                         /////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  function render(time) {
    time *= 0.02;  //timestep 20 --  convert to seconds 


    if(settings.run == 1){slime.startSimulation(slime);}



    webglUtils.resizeCanvasToDisplaySize(gl.canvas);

    gl.enable(gl.CULL_FACE);
    gl.enable(gl.DEPTH_TEST);

    // first draw from the POV of the light
    const lightWorldMatrix = m4.lookAt(
        [settings.posX, settings.posY, settings.posZ],          // position
        [settings.targetX, settings.targetY, settings.targetZ], // target
        [0, 1, 0],                                              // up
    );
    const lightProjectionMatrix = settings.perspective
        ? m4.perspective(
            degToRad(settings.fieldOfView),
            settings.projWidth / settings.projHeight,
            0.5,  // near
            100)   // far
        : m4.orthographic(
            -settings.projWidth / 2,   // left
             settings.projWidth / 2,   // right
            -settings.projHeight / 2,  // bottom
             settings.projHeight / 2,  // top
             0.5,                      // near
             100);                      // far

    // draw to the depth texture
    gl.bindFramebuffer(gl.FRAMEBUFFER, depthFramebuffer);
    gl.viewport(0, 0, depthTextureSize, depthTextureSize);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    drawScene(
        lightProjectionMatrix,
        lightWorldMatrix,
        m4.identity(),
        lightWorldMatrix,
        time,
        colorProgramInfo);

    // now draw scene to the canvas projecting the depth texture into the scene
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
    gl.clearColor(0, 0, 0, 1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    let textureMatrix = m4.identity();
    textureMatrix = m4.translate(textureMatrix, 0.5, 0.5, 0.5);
    textureMatrix = m4.scale(textureMatrix, 0.5, 0.5, 0.5);
    textureMatrix = m4.multiply(textureMatrix, lightProjectionMatrix);
    // use the inverse of this world matrix to make
    // a matrix that will transform other positions
    // to be relative this world space.
    textureMatrix = m4.multiply(
        textureMatrix,
        m4.inverse(lightWorldMatrix));

    // Compute the projection matrix
    const aspect = gl.canvas.clientWidth / gl.canvas.clientHeight;
    const projectionMatrix =
        m4.perspective(fieldOfViewRadians, aspect, 1, 2000);

    // Compute the camera's matrix using look at.
    const cameraPosition = [settings.cameraX, settings.cameraY, -(settings.cameraZ)];
    const target = [0, 0, 0];
    const up = [0, 1, 0];
    const cameraMatrix = m4.lookAt(cameraPosition, target, up);


    drawScene(
        projectionMatrix,
        cameraMatrix,
        textureMatrix,
        lightWorldMatrix,
        time,
        textureProgramInfo);


         
    // ------ Draw the frustum ------
    {
      const viewMatrix = m4.inverse(cameraMatrix);

      gl.useProgram(colorProgramInfo.program);

      // Setup all the needed attributes.
      webglUtils.setBuffersAndAttributes(gl, colorProgramInfo, cubeLinesBufferInfo);

      // scale the cube in Z so it's really long
      // to represent the texture is being projected to
      // infinity
      const mat = m4.multiply(
          lightWorldMatrix, m4.inverse(lightProjectionMatrix));

      // Set the uniforms we just computed
      webglUtils.setUniforms(colorProgramInfo, {
        u_color: [1, 1, 1, 1],
        u_view: viewMatrix,
        u_projection: projectionMatrix,
        u_world: mat,
      });

      // calls gl.drawArrays or gl.drawElements
      webglUtils.drawBufferInfo(gl, cubeLinesBufferInfo, gl.LINES);

     
      }
    requestAnimationFrame(render);
    }
  requestAnimationFrame(render);
  }

window.onload =  main();


