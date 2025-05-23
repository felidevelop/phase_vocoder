# phase_vocoder
Implementación de un Phase Vocoder clase para ajustar el tono de un sonido para un AudioWorklet

Esta es una implementacion de un phase vocoder dentro de un audio worklet llevada a cabo con el modelo OLA Processor.

Implementacion basada en el repositorio https://github.com/olvb/phaze

Toda las funciones fueron implementadas unicamente usando la clase PhaseVocoderProcessor que se extendia de AudioWorkletProcessor. Esta extension fue eliminada para manejar todas las funciones dento de la misma clase.
Unicamente la implementacion de FFT fue incluida de forma aparte.

Su uso puede verse dentro del archivo index.html y tienen un mp3 de ejemplo Totalizator.mp3, unicamente tienes que examinar un mp3, cualquiera que tengas en el pc y despues pulsar "Reproducir".

Aunque no se esta usando PHP ni Nodejs ni nada por el estilo, tienen que crear un servidor local porque la creacion de un AudioWorklet requiere de un servidor seguro https, y un localhost es considerado un servidor seguro localmente. Unicamente asi puede funcionar.

Esta implementacion tiene varias ventajas, el modelo no necesita hacer una carga previa del mp3 para ser analizado y aplicar el efecto. Ocurre en tiempo real y puedes conectar el nodo desde y hacia cualquier otro nodo de un AudioContext. Tambien puedes ajustar el pitch en cualquier momento y el sonido sera inmediatamente afectado.

Las desventajas es que este modelo a pesar de ser en tiempo real es costoso, no puedes aplicar el efecto a varios audios a la vez o puedes comsumir muchos recursos. Su implementacion en JS no es optima y por esa razon el desafio es implementar este modelo usando WASM.

De que se trata el reto?
Existe un comando llamado emcc que puede usarse en terminal. Este es capaz de compilar codigo C o C++, entre otros, y compilar el codigo a un archivo WASM. Este compilado esta a nivel de lenguaje maquina y por lo tanto su ejecucion es mucho mas optimo que usar JS.
El reto es implementar todo el modelo a C o C++, realizar una compilacion optima y crear un AudioWorklet que le encargue toda la ejecucion costoso al codigo compilado.

Logrando esta mejora optimamente, el modelo puede ejecutarse perfectamente incluso de dispositivos moviles que generalmente tienen menos potencia.
